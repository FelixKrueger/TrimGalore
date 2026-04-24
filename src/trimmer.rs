//! Trimming orchestrator for single-end and paired-end pipelines.
//!
//! This module wires together quality trimming, adapter trimming, clipping,
//! filtering, and report generation into complete processing pipelines.

use crate::alignment;
use crate::fastq::{FastqReader, FastqRecord, FastqWriter};
use crate::filters::{self, FilterResult, MaxNFilter, PairFilterResult};
use crate::quality;
use crate::report::{PairValidationStats, TrimStats};
use anyhow::Result;

/// Configuration for the trimming pipeline.
pub struct TrimConfig {
    /// R1 adapter(s): (name, sequence_bytes). Single-element for normal usage.
    pub adapters: Vec<(String, Vec<u8>)>,
    /// R2 adapter(s): (name, sequence_bytes). Empty = use R1 adapters.
    pub adapters_r2: Vec<(String, Vec<u8>)>,
    /// Max rounds of adapter trimming per read (default: 1). `-n`/`--times`.
    pub times: usize,
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

impl TrimConfig {
    /// Number of adapters used for R2 (falls back to R1 adapter count when no R2-specific adapters are set).
    pub fn r2_adapter_count(&self) -> usize {
        if self.adapters_r2.is_empty() {
            self.adapters.len()
        } else {
            self.adapters_r2.len()
        }
    }
}

/// Result of trimming a single read (carries per-read stats).
pub struct TrimResult {
    /// True if any adapter was found in any round.
    pub had_adapter: bool,
    /// Per-round matches: (adapter_index, match_len). Empty if no adapter found.
    pub adapter_matches: Vec<(usize, usize)>,
    pub quality_trimmed_bp: usize, // bases removed by quality trimming
    pub rrbs_trimmed_3prime: bool,
    pub rrbs_trimmed_5prime: bool,
    pub poly_a_trimmed: usize, // number of bases trimmed (0 = no poly-A found)
    pub poly_g_trimmed: usize, // number of bases trimmed (0 = no poly-G found)
    pub clip_5prime_applied: bool, // whether fixed 5' clip (clip_r1/clip_r2) was applied
}

/// Trim a single read in-place according to the configured pipeline.
///
/// Processing order (matches TrimGalore behavior):
/// 1. Quality trim (3' BWA algorithm)
/// 2. Adapter trim (semi-global alignment)
///    2.5. RRBS trim (MspI 2bp artifact removal, if --rrbs)
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

    // 2. Adapter trimming — multi-round best-of-N matching
    let adapters = if is_r2 && !config.adapters_r2.is_empty() {
        &config.adapters_r2
    } else {
        &config.adapters
    };

    let mut adapter_matches: Vec<(usize, usize)> = Vec::new();

    for _round in 0..config.times {
        if record.is_empty() {
            break;
        }

        let mut best: Option<(usize, alignment::AdapterMatch)> = None;
        let seq_bytes = record.seq.as_bytes();

        for (idx, (_name, adapter_seq)) in adapters.iter().enumerate() {
            if adapter_seq.is_empty() {
                continue;
            }
            if let Some(m) = alignment::find_3prime_adapter(
                seq_bytes,
                adapter_seq,
                config.error_rate,
                config.min_overlap,
            ) {
                match &best {
                    None => best = Some((idx, m)),
                    Some((_, prev)) if m.read_start < prev.read_start => {
                        best = Some((idx, m));
                    }
                    _ => {}
                }
            }
        }

        if let Some((idx, m)) = best {
            let trim_len = record.seq.len() - m.read_start;
            record.truncate(m.read_start);
            adapter_matches.push((idx, trim_len));
        } else {
            break;
        }
    }

    let had_adapter = !adapter_matches.is_empty();

    // 2.5. RRBS trimming (MspI 2bp end-repair artifact removal)
    let mut rrbs_trimmed_3prime = false;
    let mut rrbs_trimmed_5prime = false;

    if config.rrbs {
        if config.non_directional {
            // Non-directional: check all reads (R1 and R2)
            if record.seq.len() > 2 {
                let seq_bytes = record.seq.as_bytes();
                if seq_bytes.len() >= 3 && (seq_bytes[..3] == *b"CAA" || seq_bytes[..3] == *b"CGA")
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
            if !(config.is_paired && is_r2) && record.seq.len() >= 2 && had_adapter {
                record.truncate(record.seq.len() - 2);
                rrbs_trimmed_3prime = true;
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

    // 3. Trim Ns from both ends.
    // Suppressed under --rrbs to match Perl v0.6.x: Perl's RRBS code path
    // omits `$trim_n` from its Cutadapt invocations (trim_galore:876–915,
    // 1353/1358/1365 in the non-RRBS path for contrast). Keeps v2.x output
    // byte-identical to v0.6.x for --trim-n + --rrbs users.
    if config.trim_n && !config.rrbs {
        record.trim_ns();
    }

    // 4. Fixed clipping
    let clip_5 = if is_r2 {
        config.clip_r2
    } else {
        config.clip_r1
    };
    let clip_3 = if is_r2 {
        config.three_prime_clip_r2
    } else {
        config.three_prime_clip_r1
    };

    let mut clip_5prime_applied = false;
    if let Some(n) = clip_5 {
        let clipped = record.clip_5prime(n);
        clip_5prime_applied = clipped.is_some();
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
        adapter_matches,
        quality_trimmed_bp,
        rrbs_trimmed_3prime,
        rrbs_trimmed_5prime,
        poly_a_trimmed,
        poly_g_trimmed,
        clip_5prime_applied,
    }
}

/// Update per-adapter stats from a trim result (handles multi-round matches).
pub fn update_adapter_stats(stats: &mut TrimStats, result: &TrimResult) {
    if result.had_adapter {
        stats.total_reads_with_adapter += 1;
        for &(idx, len) in &result.adapter_matches {
            stats.reads_with_adapter[idx] += 1;
            if len >= stats.adapter_length_counts[idx].len() {
                stats.adapter_length_counts[idx].resize(len + 1, 0);
            }
            stats.adapter_length_counts[idx][len] += 1;
        }
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
    let mut stats = TrimStats::with_adapter_count(config.adapters.len());

    while let Some(mut record) = reader.next_record()? {
        stats.total_reads += 1;
        stats.total_bp_processed += record.seq.len();

        let result = trim_read(&mut record, config, false);
        stats.bases_quality_trimmed += result.quality_trimmed_bp;
        update_adapter_stats(&mut stats, &result);
        if result.rrbs_trimmed_3prime {
            stats.rrbs_trimmed_3prime += 1;
        }
        if result.rrbs_trimmed_5prime {
            stats.rrbs_trimmed_5prime += 1;
        }
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

        // Track bp after trimming but before length filtering (Cutadapt-compatible stat)
        stats.total_bp_after_trim += record.seq.len();

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
#[allow(clippy::too_many_arguments)]
pub fn run_paired_end(
    reader_r1: &mut FastqReader,
    reader_r2: &mut FastqReader,
    writer_r1: &mut FastqWriter,
    writer_r2: &mut FastqWriter,
    mut unpaired_r1: Option<&mut FastqWriter>,
    mut unpaired_r2: Option<&mut FastqWriter>,
    config: &TrimConfig,
    unpaired: crate::filters::UnpairedLengths,
) -> Result<(TrimStats, TrimStats, PairValidationStats)> {
    let mut stats_r1 = TrimStats::with_adapter_count(config.adapters.len());
    let mut stats_r2 = TrimStats::with_adapter_count(config.r2_adapter_count());
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
                update_adapter_stats(&mut stats_r1, &result_r1);
                update_adapter_stats(&mut stats_r2, &result_r2);
                if result_r1.rrbs_trimmed_3prime {
                    stats_r1.rrbs_trimmed_3prime += 1;
                }
                if result_r1.rrbs_trimmed_5prime {
                    stats_r1.rrbs_trimmed_5prime += 1;
                }
                if result_r2.rrbs_trimmed_3prime {
                    stats_r2.rrbs_trimmed_3prime += 1;
                }
                if result_r2.rrbs_trimmed_5prime {
                    stats_r2.rrbs_trimmed_5prime += 1;
                }
                if config.rrbs && result_r2.clip_5prime_applied {
                    stats_r2.rrbs_r2_clipped_5prime += 1;
                }
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

                // Track bp after trimming but before pair/length filtering (Cutadapt-compatible stat)
                stats_r1.total_bp_after_trim += r1.seq.len();
                stats_r2.total_bp_after_trim += r2.seq.len();

                // Pair-aware filtering
                match filters::filter_paired_end(
                    &r1,
                    &r2,
                    config.length_cutoff,
                    config.max_length,
                    config.max_n.clone(),
                    unpaired,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::FastqRecord;

    /// Build a minimal TrimConfig for adapter-trimming tests.
    /// Quality cutoff 0 disables quality trimming so we test only adapter logic.
    fn test_config_with_adapters(adapters: Vec<(&str, &str)>, times: usize) -> TrimConfig {
        TrimConfig {
            adapters: adapters
                .iter()
                .map(|(name, seq)| (name.to_string(), seq.as_bytes().to_vec()))
                .collect(),
            adapters_r2: Vec::new(),
            times,
            quality_cutoff: 0,
            phred_offset: 33,
            error_rate: 0.1,
            min_overlap: 1,
            length_cutoff: 0,
            max_length: None,
            max_n: None,
            trim_n: false,
            clip_r1: None,
            clip_r2: None,
            three_prime_clip_r1: None,
            three_prime_clip_r2: None,
            rename: false,
            nextseq: false,
            rrbs: false,
            non_directional: false,
            is_paired: false,
            poly_a: false,
            poly_g: false,
            discard_untrimmed: false,
        }
    }

    fn make_record(seq: &str) -> FastqRecord {
        FastqRecord {
            id: "@test_read".to_string(),
            seq: seq.to_string(),
            qual: "I".repeat(seq.len()),
        }
    }

    // ── Multi-adapter best-of-N tests ───────────────────────────

    #[test]
    fn test_trim_read_multi_adapter_best_wins() {
        // Read: 20bp random + CTGTCTCTTATA (Nextera, 12bp) + 8bp padding + AGATCGGAAGAGC (Illumina, 13bp)
        // Both adapters are present. Nextera starts earlier (position 20) so trims more → it should win.
        let read_seq = "ACGTACGTACGTACGTACGTCTGTCTCTTATAGGGGGGGGAGATCGGAAGAGC";
        //              |----20bp random----|--Nextera 12bp--|--8bp pad--|--Illumina 13bp--|
        let config = test_config_with_adapters(
            vec![("illumina", "AGATCGGAAGAGC"), ("nextera", "CTGTCTCTTATA")],
            1,
        );

        let mut record = make_record(read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        assert_eq!(result.adapter_matches.len(), 1);
        // Nextera (index 1) should win because it starts at position 20 (trims more)
        assert_eq!(result.adapter_matches[0].0, 1); // adapter index = nextera
        assert_eq!(record.seq.len(), 20); // trimmed at position 20
    }

    #[test]
    fn test_trim_read_multi_adapter_single_match() {
        // Read has only the Illumina adapter, not Nextera.
        let read_seq = "ACGTACGTACGTACGTACGTACGTACGTAAGAGATCGGAAGAGC";
        //              |------30bp random + partial-------|--Illumina--|
        let config = test_config_with_adapters(
            vec![
                ("illumina", "AGATCGGAAGAGC"),
                ("nextera", "CTGTCTCTTATA"),
                ("smallrna", "TGGAATTCTCGG"),
            ],
            1,
        );

        let mut record = make_record(read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        assert_eq!(result.adapter_matches.len(), 1);
        assert_eq!(result.adapter_matches[0].0, 0); // illumina is index 0
    }

    #[test]
    fn test_trim_read_multi_adapter_none_match() {
        // Read with no adapter sequences at all.
        let read_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let config = test_config_with_adapters(
            vec![("illumina", "AGATCGGAAGAGC"), ("nextera", "CTGTCTCTTATA")],
            1,
        );

        let mut record = make_record(read_seq);
        let original_len = record.seq.len();
        let result = trim_read(&mut record, &config, false);

        assert!(!result.had_adapter);
        assert!(result.adapter_matches.is_empty());
        assert_eq!(record.seq.len(), original_len); // unchanged
    }

    #[test]
    fn test_trim_read_single_adapter_unchanged() {
        // Regression test: single adapter behavior identical to before multi-adapter.
        let read_seq = "ACGTACGTACGTACGTACGTAAGAGATCGGAAGAGC";
        let config = test_config_with_adapters(vec![("illumina", "AGATCGGAAGAGC")], 1);

        let mut record = make_record(read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        assert_eq!(result.adapter_matches.len(), 1);
        assert_eq!(result.adapter_matches[0].0, 0); // adapter index 0
        // The adapter starts matching somewhere in the read; verify the read got shorter
        assert!(record.seq.len() < read_seq.len());
    }

    // ── Multi-round (-n / --times) tests ────────────────────────

    #[test]
    fn test_trim_read_times_2() {
        // Read: 20bp random + Nextera adapter + Illumina adapter (concatenated at 3' end).
        // With -n 2 and BOTH adapters specified:
        //   Round 1: best-of-N picks Nextera (starts at pos 20, trims more than Illumina at pos 32)
        //            → truncates at pos 20, remaining = "ACGTACGTACGTACGTACGT"
        //   Round 2: remaining 20bp has no adapter → stops.
        //
        // Alternative scenario with single adapter and -n 2:
        // Use the same adapter twice in the read to test multi-round.
        let adapter = "AGATCGGAAGAGC";
        let read_seq = format!("ACGTACGTACGTACGT{}{}", adapter, adapter);
        // = "ACGTACGTACGTACGTAGATCGGAAGAGCAGATCGGAAGAGC" (16 + 13 + 13 = 42bp)
        let config = test_config_with_adapters(vec![("illumina", adapter)], 2);

        let mut record = make_record(&read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        // Round 1: adapter found at position 16 → trims to "ACGTACGTACGTACGT" (16bp)
        // Round 2: 16bp with no adapter → stops
        // So only 1 match (the first round trims BOTH copies at once since
        // the alignment finds the leftmost match).
        assert!(!result.adapter_matches.is_empty());
        assert_eq!(record.seq, "ACGTACGTACGTACGT");
    }

    #[test]
    fn test_trim_read_times_stops_on_no_match() {
        // -n 3, but only one adapter present → should stop after 1 productive round.
        // Use min_overlap=3 (realistic default) and a read whose remaining portion
        // after trimming is pure TTTT — no chance of accidental adapter matches.
        let read_seq = "TTTTTTTTTTTTTTTTTTTTAAGAGATCGGAAGAGC";
        let mut config = test_config_with_adapters(vec![("illumina", "AGATCGGAAGAGC")], 3);
        config.min_overlap = 3;

        let mut record = make_record(read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        // Only 1 round should match — remaining TTTs can't match the adapter
        assert_eq!(result.adapter_matches.len(), 1);
        // Verify the read was trimmed to just the T stretch
        assert!(!record.seq.contains("AGATCG"));
    }

    #[test]
    fn test_trim_read_times_1_default() {
        // -n 1 (default): even if the read has adapter appearing twice,
        // only one round of trimming occurs.
        // Construct: random + adapter + adapter (concatenated adapters)
        let adapter = "AGATCGGAAGAGC";
        let read_seq = format!("ACGTACGTACGTACGT{}{}", adapter, adapter);
        let config = test_config_with_adapters(
            vec![("illumina", adapter)],
            1, // default: single round
        );

        let mut record = make_record(&read_seq);
        let result = trim_read(&mut record, &config, false);

        assert!(result.had_adapter);
        // Only 1 match despite adapter appearing twice
        assert_eq!(result.adapter_matches.len(), 1);
        // The first (leftmost) match wins, trimming at position 16
        assert_eq!(record.seq, "ACGTACGTACGTACGT");
    }

    // ── --trim-n × --rrbs interaction (Perl v0.6.x parity) ──────

    #[test]
    fn test_trim_n_without_rrbs_trims_trailing_ns() {
        let mut config = test_config_with_adapters(vec![], 1);
        config.trim_n = true;
        config.rrbs = false;

        let mut record = make_record("ACGTACGTACGTNNNN");
        trim_read(&mut record, &config, false);
        assert_eq!(record.seq, "ACGTACGTACGT"); // 4 trailing Ns trimmed
    }

    #[test]
    fn test_trim_n_with_rrbs_suppressed() {
        // Perl v0.6.x parity: --trim-n is a no-op under --rrbs because Perl's
        // RRBS code path omits $trim_n from its Cutadapt invocations.
        let mut config = test_config_with_adapters(vec![], 1);
        config.trim_n = true;
        config.rrbs = true;

        let mut record = make_record("ACGTACGTACGTNNNN");
        trim_read(&mut record, &config, false);
        assert_eq!(record.seq, "ACGTACGTACGTNNNN"); // Ns preserved
    }

    #[test]
    fn test_trim_n_rrbs_also_preserves_leading_ns() {
        // trim_ns() handles both ends — verify the RRBS suppression covers
        // the leading-N case too, not just trailing.
        let mut config = test_config_with_adapters(vec![], 1);
        config.trim_n = true;
        config.rrbs = true;

        let mut record = make_record("NNACGTACGT");
        trim_read(&mut record, &config, false);
        assert_eq!(record.seq, "NNACGTACGT");
    }
}
