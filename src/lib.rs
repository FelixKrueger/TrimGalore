//! Trim Galore's library target exists so the binary and the integration tests can
//! reach internals. Its surface is an implementation detail and carries no semver
//! guarantee: modules are `pub` for reachability, not as a published API (#402).

pub mod adapter;
pub mod alignment;
pub mod bam;
pub mod cli;
pub mod clump;
pub mod clump_only;
pub mod demux;
pub mod fastq;
pub mod fastqc;
pub mod filters;
pub mod format;
pub mod io;
pub mod library;
pub mod parallel;
pub mod quality;
pub mod report;
pub mod specialty;
pub mod trimmer;
