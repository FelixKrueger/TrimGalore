//! Trimming orchestrator for single-end and paired-end pipelines.
//!
//! This module wires together quality trimming, adapter trimming, clipping,
//! filtering, and report generation into complete processing pipelines.

use anyhow::Result;
use crate::alignment;
use crate::fastq::{FastqReader, FastqRecord, FastqWriter};
use crate::filters::{self, FilterResult, MaxNFilter, PairFilterResult};
use crate::quality;
use crate::report::{TrimStats, PairValidationStats};

/// Configuration for the trimming pipeline.
pub struct TrimConfig {
    pub adapter: Vec<u8>,
    pub adapter_r2: Option<Vec<u8>>,
    pub quality_cutoff: u8,
    pub phred_offset: u8,
    pub error_rate: f64,
    pub min_overlap: usize,
    pub length_cutoff: usize,
    pub max_length: Option<usize>,
    pub max_n: Option<MaxNFilter>,
    pub trim_n: bool,
    pub clip_r1: Option<usize>,
    pub clip_r2: Option<usize>,
    pub three_prime_clip_r1: Option<usize>,
    pub three_prime_clip_r2: Option<usize>,
    pub rename: bool,
    pub nextseq: bool,
    pub rrbs: bool,
    pub non_directional: bool,
    pub is_paired: bool,
    pub poly_a: bool,
    pub poly_g: bool,
    pub discard_untrimmed: bool,
}

/// Result of trimming a single read (carries per-read stats).
pub struct TrimResult {
    pub had_adapter: bool,
    pub adapter_match_len: usize, // length of adapter overlap found (0 = none)
    pub quality_trimmed_bp: usize, // bases removed by quality trimming
    pub rrbs_trimmed_3prime: bool,
    pub rrbs_trimmed_5prime: bool,
    pub poly_a_trimmed: usize, // number of bases trimmed (0 = no poly-A found)
    pub poly_g_trimmed: usize, // number of bases trimmed (0 = no poly-G found)
}

/// Trim a single read in-place according to the configured pipeline.
///
/// Processing order (matches TrimGalore behavior):
/// 1. Quality trim (3' BWA algorithm)
/// 2. Adapter trim (semi-global alignment)
/// 2.5. RRBS trim (MspI 2bp artifact removal, if --rrbs)
/// 3. Trim Ns (if --trim-n)
/// 4. 5'/3' fixed clipping
///
/// Returns a TrimResult with adapter and RRBS status.
pub fn trim_read(record: &mut FastqRecord, config: &TrimConfig, is_r2: bool) -> TrimResult {
    // 1. Quality trimming (NextSeq mode overrides G-base quality to 0)
    let len_before_quality = record.seq.len();
    if config.quality_cutoff > 0 {
        let trim_pos = if config.nextseq {
            quality::quality_trim_3prime_nextseq(
                record.seq.as_bytes(),
                record.qual.as_bytes(),
                config.quality_cutoff,
                config.phred_offset,
            )
        } else {
            quality::quality_trim_3prime(
                record.qual.as_bytes(),
                config.quality_cutoff,
                config.phred_offset,
            )
        };
        record.truncate(trim_pos);
    }
    let quality_trimmed_bp = len_before_quality - record.seq.len();

    // 2. Adapter trimming
    let adapter = if is_r2 {
        config.adapter_r2.as_deref().unwrap_or(&config.adapter)
    } else {
        &config.adapter
    };

    let mut adapter_match_len: usize = 0;
    let had_adapter = if !adapter.is_empty() && !record.is_empty() {
        let seq_len = record.seq.len();
        if let Some(m) = alignment::find_3prime_adapter(
            record.seq.as_bytes(),
            adapter,
            config.error_rate,
            config.min_overlap,
        ) {
            adapter_match_len = seq_len - m.read_start;
            record.truncate(m.read_start);
            true
        } else {
            false
        }
    } else {
        false
    };

    // 2.5. RRBS trimming (MspI 2bp end-repair artifact removal)
    let mut rrbs_trimmed_3prime = false;
    let mut rrbs_trimmed_5prime = false;

    if config.rrbs {
        if config.non_directional {
            // Non-directional: check all reads (R1 and R2)
            if record.seq.len() > 2 {
                let seq_bytes = record.seq.as_bytes();
                if seq_bytes.len() >= 3
                    && (seq_bytes[..3] == *b"CAA" || seq_bytes[..3] == *b"CGA")
                {
                    // Non-directional artifact: trim 2bp from 5' end
                    record.clip_5prime(2);
                    rrbs_trimmed_5prime = true;
                } else if had_adapter {
                    // Standard directional-style 3' trim
                    record.truncate(record.seq.len() - 2);
                    rrbs_trimmed_3prime = true;
                }
            }
        } else {
            // Directional RRBS: only R1/SE get 3' trimmed; R2 handled by auto-set clip_r2=2
            if !(config.is_paired && is_r2) {
                if record.seq.len() >= 2 && had_adapter {
                    record.truncate(record.seq.len() - 2);
                    rrbs_trimmed_3prime = true;
                }
            }
        }
    }

    // 2.7. Poly-A / Poly-T trimming (after adapter + RRBS, before N-trimming)
    // R1/SE: trim poly-A from 3' end. R2: trim poly-T from 5' end.
    let mut poly_a_trimmed: usize = 0;
    if config.poly_a && !record.is_empty() {
        let revcomp = is_r2; // R2 gets poly-T (5' end) trimming
        let seq_len_before = record.seq.len();
        let idx = quality::poly_a_trim_index(record.seq.as_bytes(), revcomp);
        if revcomp {
            // Poly-T head: clip idx bases from 5' end
            if idx > 0 {
                record.clip_5prime(idx);
                poly_a_trimmed = idx;
            }
        } else {
            // Poly-A tail: truncate at idx
            if idx < seq_len_before {
                record.truncate(idx);
                poly_a_trimmed = seq_len_before - idx;
            }
        }
    }

    // 2.8. Poly-G / Poly-C trimming (2-colour chemistry artifact removal)
    // R1/SE: trim poly-G from 3' end. R2: trim poly-C from 5' end.
    let mut poly_g_trimmed: usize = 0;
    if config.poly_g && !record.is_empty() {
        let revcomp = is_r2;
        let seq_len_before = record.seq.len();
        let idx = quality::homopolymer_trim_index(record.seq.as_bytes(), b'G', revcomp);
        if revcomp {
            if idx > 0 {
                record.clip_5prime(idx);
                poly_g_trimmed = idx;
            }
        } else {
            if idx < seq_len_before {
                record.truncate(idx);
                poly_g_trimmed = seq_len_before - idx;
            }
        }
    }

    // 3. Trim Ns from both ends
    if config.trim_n {
        record.trim_ns();
    }

    // 4. Fixed clipping
    let clip_5 = if is_r2 { config.clip_r2 } else { config.clip_r1 };
    let clip_3 = if is_r2 { config.three_prime_clip_r2 } else { config.three_prime_clip_r1 };

    if let Some(n) = clip_5 {
        let clipped = record.clip_5prime(n);
        if config.rename {
            if let Some(seq) = clipped {
                record.append_to_id(&format!(":clip5:{}", seq));
            }
        }
    }

    if let Some(n) = clip_3 {
        let clipped = record.clip_3prime(n);
        if config.rename {
            if let Some(seq) = clipped {
                record.append_to_id(&format!(":clip3:{}", seq));
            }
        }
    }

    TrimResult {
        had_adapter,
        adapter_match_len,
        quality_trimmed_bp,
        rrbs_trimmed_3prime,
        rrbs_trimmed_5prime,
        poly_a_trimmed,
        poly_g_trimmed,
    }
}

/// Run the single-end trimming pipeline.
///
/// Reads input, trims each read, applies filters, writes output.
/// Returns trimming statistics.
pub fn run_single_end(
    reader: &mut FastqReader,
    writer: &mut FastqWriter,
    config: &TrimConfig,
) -> Result<TrimStats> {
    let mut stats = TrimStats::default();

    while let Some(mut record) = reader.next_record()? {
        stats.total_reads += 1;
        stats.total_bp_processed += record.seq.len();

        let result = trim_read(&mut record, config, false);
        stats.bases_quality_trimmed += result.quality_trimmed_bp;
        if result.had_adapter {
            stats.reads_with_adapter += 1;
            let len = result.adapter_match_len;
            if len >= stats.adapter_length_counts.len() {
                stats.adapter_length_counts.resize(len + 1, 0);
            }
            stats.adapter_length_counts[len] += 1;
        }
        if result.rrbs_trimmed_3prime { stats.rrbs_trimmed_3prime += 1; }
        if result.rrbs_trimmed_5prime { stats.rrbs_trimmed_5prime += 1; }
        if result.poly_a_trimmed > 0 {
            stats.poly_a_trimmed += 1;
            stats.poly_a_bases_trimmed += result.poly_a_trimmed;
        }
        if result.poly_g_trimmed > 0 {
            stats.poly_g_trimmed += 1;
            stats.poly_g_bases_trimmed += result.poly_g_trimmed;
        }

        // Discard reads without adapter match (--discard-untrimmed)
        if config.discard_untrimmed && !result.had_adapter {
            stats.discarded_untrimmed += 1;
            continue;
        }

        // Apply filters
        match filters::filter_single_end(
            &record,
            config.length_cutoff,
            config.max_length,
            config.max_n.clone(),
        ) {
            FilterResult::Pass => {
                stats.total_bp_written += record.seq.len();
                writer.write_record(&record)?;
                stats.reads_written += 1;
            }
            FilterResult::TooShort => stats.too_short += 1,
            FilterResult::TooLong => stats.too_long += 1,
            FilterResult::TooManyN => stats.too_many_n += 1,
        }
    }

    Ok(stats)
}

/// Run the paired-end trimming pipeline (single-pass).
///
/// Reads R1 and R2 in lockstep, trims both, applies pair-aware filtering.
/// This is the key architectural improvement over TrimGalore's 3-pass approach.
pub fn run_paired_end(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    writer_r1: &mut FastqWriter,
    writer_r2: &mut FastqWriter,
    mut unpaired_r1: Option<&mut FastqWriter>,
    mut unpaired_r2: Option<&mut FastqWriter>,
    config: &TrimConfig,
    unpaired_length_r1: usize,
    unpaired_length_r2: usize,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    let mut stats_r1 = TrimStats::default();
    let mut stats_r2 = TrimStats::default();
    let mut pair_stats = PairValidationStats::default();

    loop {
        let rec1 = reader_r1.next_record()?;
        let rec2 = reader_r2.next_record()?;

        match (rec1, rec2) {
            (Some(mut r1), Some(mut r2)) => {
                stats_r1.total_reads += 1;
                stats_r2.total_reads += 1;
                stats_r1.total_bp_processed += r1.seq.len();
                stats_r2.total_bp_processed += r2.seq.len();
                pair_stats.pairs_analyzed += 1;

                // Trim both reads
                let result_r1 = trim_read(&mut r1, config, false);
                let result_r2 = trim_read(&mut r2, config, true);

                stats_r1.bases_quality_trimmed += result_r1.quality_trimmed_bp;
                stats_r2.bases_quality_trimmed += result_r2.quality_trimmed_bp;
                if result_r1.had_adapter {
                    stats_r1.reads_with_adapter += 1;
                    let len = result_r1.adapter_match_len;
                    if len >= stats_r1.adapter_length_counts.len() {
                        stats_r1.adapter_length_counts.resize(len + 1, 0);
                    }
                    stats_r1.adapter_length_counts[len] += 1;
                }
                if result_r2.had_adapter {
                    stats_r2.reads_with_adapter += 1;
                    let len = result_r2.adapter_match_len;
                    if len >= stats_r2.adapter_length_counts.len() {
                        stats_r2.adapter_length_counts.resize(len + 1, 0);
                    }
                    stats_r2.adapter_length_counts[len] += 1;
                }
                if result_r1.rrbs_trimmed_3prime { stats_r1.rrbs_trimmed_3prime += 1; }
                if result_r1.rrbs_trimmed_5prime { stats_r1.rrbs_trimmed_5prime += 1; }
                if result_r2.rrbs_trimmed_3prime { stats_r2.rrbs_trimmed_3prime += 1; }
                if result_r2.rrbs_trimmed_5prime { stats_r2.rrbs_trimmed_5prime += 1; }
                if result_r1.poly_a_trimmed > 0 {
                    stats_r1.poly_a_trimmed += 1;
                    stats_r1.poly_a_bases_trimmed += result_r1.poly_a_trimmed;
                }
                if result_r2.poly_a_trimmed > 0 {
                    stats_r2.poly_a_trimmed += 1;
                    stats_r2.poly_a_bases_trimmed += result_r2.poly_a_trimmed;
                }
                if result_r1.poly_g_trimmed > 0 {
                    stats_r1.poly_g_trimmed += 1;
                    stats_r1.poly_g_bases_trimmed += result_r1.poly_g_trimmed;
                }
                if result_r2.poly_g_trimmed > 0 {
                    stats_r2.poly_g_trimmed += 1;
                    stats_r2.poly_g_bases_trimmed += result_r2.poly_g_trimmed;
                }

                // Discard pairs without adapter match (--discard-untrimmed)
                if config.discard_untrimmed && !result_r1.had_adapter && !result_r2.had_adapter {
                    stats_r1.discarded_untrimmed += 1;
                    stats_r2.discarded_untrimmed += 1;
                    pair_stats.pairs_removed += 1;
                    continue;
                }

                // Pair-aware filtering
                match filters::filter_paired_end(
                    &r1,
                    &r2,
                    config.length_cutoff,
                    config.max_length,
                    config.max_n.clone(),
                    unpaired_length_r1,
                    unpaired_length_r2,
                ) {
                    PairFilterResult::Pass => {
                        stats_r1.total_bp_written += r1.seq.len();
                        stats_r2.total_bp_written += r2.seq.len();
                        writer_r1.write_record(&r1)?;
                        writer_r2.write_record(&r2)?;
                        stats_r1.reads_written += 1;
                        stats_r2.reads_written += 1;
                    }
                    PairFilterResult::TooManyN => {
                        // N-filter: always discard entire pair, no rescue
                        pair_stats.pairs_removed += 1;
                        pair_stats.pairs_removed_n += 1;
                        stats_r1.too_many_n += 1;
                        stats_r2.too_many_n += 1;
                    }
                    PairFilterResult::TooShort { r1_ok, r2_ok } => {
                        pair_stats.pairs_removed += 1;
                        stats_r1.too_short += 1;
                        stats_r2.too_short += 1;

                        // Rescue individual reads if --retain_unpaired
                        if let Some(ref mut w) = unpaired_r1.as_deref_mut() {
                            if r1_ok {
                                w.write_record(&r1)?;
                                pair_stats.r1_unpaired += 1;
                            }
                        }
                        if let Some(ref mut w) = unpaired_r2.as_deref_mut() {
                            if r2_ok {
                                w.write_record(&r2)?;
                                pair_stats.r2_unpaired += 1;
                            }
                        }
                    }
                    PairFilterResult::TooLong => {
                        pair_stats.pairs_removed += 1;
                        pair_stats.pairs_removed_too_long += 1;
                        stats_r1.too_long += 1;
                        stats_r2.too_long += 1;
                    }
                }
            }
            (None, None) => break, // Both files ended
            (Some(_), None) => {
                anyhow::bail!(
                    "Read 2 file is truncated — R1 has more reads than R2. \
                     Please check your paired-end input files!"
                );
            }
            (None, Some(_)) => {
                anyhow::bail!(
                    "Read 1 file is truncated — R2 has more reads than R1. \
                     Please check your paired-end input files!"
                );
            }
        }
    }

    Ok((stats_r1, stats_r2, pair_stats))
}
