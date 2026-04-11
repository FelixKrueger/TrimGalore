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
}

/// Trim a single read in-place according to the configured pipeline.
///
/// Processing order (matches TrimGalore behavior):
/// 1. Quality trim (3' BWA algorithm)
/// 2. Adapter trim (semi-global alignment)
/// 3. Trim Ns (if --trim-n)
/// 4. 5'/3' fixed clipping
///
/// Returns whether an adapter was found.
pub fn trim_read(record: &mut FastqRecord, config: &TrimConfig, is_r2: bool) -> bool {
    // 1. Quality trimming
    if config.quality_cutoff > 0 {
        let trim_pos = quality::quality_trim_3prime(
            record.qual.as_bytes(),
            config.quality_cutoff,
            config.phred_offset,
        );
        record.truncate(trim_pos);
    }

    // 2. Adapter trimming
    let adapter = if is_r2 {
        config.adapter_r2.as_deref().unwrap_or(&config.adapter)
    } else {
        &config.adapter
    };

    let had_adapter = if !adapter.is_empty() && !record.is_empty() {
        if let Some(m) = alignment::find_3prime_adapter(
            record.seq.as_bytes(),
            adapter,
            config.error_rate,
            config.min_overlap,
        ) {
            record.truncate(m.read_start);
            true
        } else {
            false
        }
    } else {
        false
    };

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

    had_adapter
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

        let had_adapter = trim_read(&mut record, config, false);
        if had_adapter {
            stats.reads_with_adapter += 1;
        }

        // Apply filters
        match filters::filter_single_end(
            &record,
            config.length_cutoff,
            config.max_length,
            config.max_n.clone(),
        ) {
            FilterResult::Pass => {
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
                pair_stats.pairs_analyzed += 1;

                // Trim both reads
                let adapter_r1 = trim_read(&mut r1, config, false);
                let adapter_r2 = trim_read(&mut r2, config, true);

                if adapter_r1 { stats_r1.reads_with_adapter += 1; }
                if adapter_r2 { stats_r2.reads_with_adapter += 1; }

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
