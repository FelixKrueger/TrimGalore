//! Reference sketch — paired-uBAM streaming de-interleaver with bounded-buffer
//! adversarial-input detection.
//!
//! NOT compiled or wired into the crate. Carry the algorithm into `src/bam.rs`
//! when implementing the plan.
//!
//! Strategy:
//!   * Strictly-interleaved happy path (`R1, R2, R1, R2, …`) → drains in O(1) memory.
//!   * Near-interleaved slack (a small number of consecutive same-side records,
//!     e.g. one R1 arriving before its R2 because of a sort tie-break) → tolerated
//!     up to MAX_SLACK pending records on either side.
//!   * Truly-grouped input (`R1×N, R2×N`) → as soon as the pending side queue
//!     exceeds MAX_SLACK, error out with a clear remediation message. Bounded
//!     memory: MAX_SLACK × sizeof(FastqRecord), not the whole file.
//!
//! Recommended MAX_SLACK = 1024 (~ 1 MB at 1 KB/record). See SPIKE report.

use std::collections::VecDeque;

const MAX_SLACK: usize = 1024;

/// Outcome of attempting to feed one BAM record into the de-interleaver.
pub enum Feed {
    /// A matched pair `(r1, r2)` is now available. Push to consumers.
    Pair { r1: FastqRecord, r2: FastqRecord },
    /// Record buffered; no pair yet.
    Buffered,
    /// Either side queue exceeded MAX_SLACK — the input is grouped, not
    /// interleaved. Caller propagates this as a fatal error.
    GroupedInputDetected { pending_r1: usize, pending_r2: usize },
    /// Flag is neither R1 nor R2 (or both), or zero. Caller errors.
    InvalidFlag { flag: u16 },
}

pub struct DeInterleaver {
    r1_queue: VecDeque<FastqRecord>,
    r2_queue: VecDeque<FastqRecord>,
}

impl DeInterleaver {
    pub fn new() -> Self {
        Self {
            r1_queue: VecDeque::new(),
            r2_queue: VecDeque::new(),
        }
    }

    /// Feed one converted BAM record (already `FastqRecord` + side flag bits).
    pub fn feed(&mut self, rec: FastqRecord, flags: u16) -> Feed {
        const FREAD1: u16 = 0x40;
        const FREAD2: u16 = 0x80;
        let is_r1 = flags & FREAD1 != 0;
        let is_r2 = flags & FREAD2 != 0;

        match (is_r1, is_r2) {
            (true, false) => self.r1_queue.push_back(rec),
            (false, true) => self.r2_queue.push_back(rec),
            _ => return Feed::InvalidFlag { flag: flags },
        }

        // Emit a pair if both queues are non-empty (mate-name pairing is
        // checked separately — see name_check below; in strict-interleaving
        // mode they will pop in order).
        if !self.r1_queue.is_empty() && !self.r2_queue.is_empty() {
            let r1 = self.r1_queue.pop_front().unwrap();
            let r2 = self.r2_queue.pop_front().unwrap();
            // Optional sanity: r1.id == r2.id (after /1 /2 stripping)
            return Feed::Pair { r1, r2 };
        }

        if self.r1_queue.len() > MAX_SLACK || self.r2_queue.len() > MAX_SLACK {
            return Feed::GroupedInputDetected {
                pending_r1: self.r1_queue.len(),
                pending_r2: self.r2_queue.len(),
            };
        }

        Feed::Buffered
    }

    /// Call at EOF; returns Err if any unmatched records remain.
    pub fn finish(self) -> anyhow::Result<()> {
        if !self.r1_queue.is_empty() || !self.r2_queue.is_empty() {
            anyhow::bail!(
                "paired uBAM ended with {} unmatched R1 and {} unmatched R2 records",
                self.r1_queue.len(),
                self.r2_queue.len()
            );
        }
        Ok(())
    }
}

/// Suggested error text when `GroupedInputDetected` fires:
const GROUPED_INPUT_ERR: &str =
    "Paired uBAM input is grouped (all R1 records arrive before R2 records).\n\
     TrimGalore's paired-uBAM reader requires mate-adjacent ordering to stream\n\
     without unbounded buffering.\n\
     \n\
     Re-interleave the input first, then re-run TrimGalore:\n\
     \n\
         samtools collate -O input.bam tmp_prefix > interleaved.bam\n\
     \n\
     (or `samtools sort -n` if you prefer a fully name-sorted file).";

// Placeholder type — real code uses crate::fastq::FastqRecord.
pub struct FastqRecord;
