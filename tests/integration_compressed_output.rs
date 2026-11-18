//! The **compressed** bytes of gzip output, pinned.
//!
//! Everything else here compares output through `gzip -dc`, and the CI validation
//! job does the same — so nothing could see a change in the gzip framing itself.
//! Downstream can: nf-core/modules pins the compressed md5s in its nf-test
//! snapshot, so a change in the teardown's flushes or the compression level breaks
//! a consumer without failing a single test of ours.
//!
//! `fastq::PRE_434_GZ_SYNC_FLUSHES` exists for that contract, and the test that
//! used to guard it inferred the flush count from the bytes — which stopped
//! working when flate2 1.1.11 began coalescing sync flushes across an output
//! buffer boundary, making one flush and two produce identical output. These
//! assertions infer nothing; they pin the bytes a consumer sees.
//!
//! **If one of these fails, do not just re-pin it.** Find what changed the output
//! and whether downstream can take it. A deliberate change means re-pinning here
//! *and* telling nf-core/modules.

use std::path::Path;
use std::process::Command;

use tempfile::tempdir;

/// CRC32 and length of a file's bytes. Reproduce by hand with:
/// `python3 -c "import zlib,sys; b=open(sys.argv[1],'rb').read(); print(hex(zlib.crc32(b)), len(b))" FILE`
fn digest(path: &Path) -> (u32, usize) {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    let mut h = crc32fast::Hasher::new();
    h.update(&bytes);
    (h.finalize(), bytes.len())
}

fn trim(dir: &Path, extra: &[&str]) -> std::path::PathBuf {
    let mut args: Vec<&str> = vec!["--output_dir", dir.to_str().unwrap()];
    args.extend_from_slice(extra);
    args.push("test_files/illumina_10K.fastq.gz");
    let out = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args(&args)
        .output()
        .expect("trim_galore runs");
    assert!(
        out.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    dir.join("illumina_10K_trimmed.fq.gz")
}

/// Both values verified identical under flate2 1.1.10 and 1.1.11, which differ in
/// how they coalesce the teardown's sync flushes — so these pin Trim Galore's
/// output rather than the library's internals.
const SERIAL_GZ: (u32, usize) = (0xc8e4_198e, 1_447_095);
const PARALLEL_GZ: (u32, usize) = (0x9def_f60b, 1_445_466);

#[test]
fn serial_gzip_output_bytes_are_pinned() {
    let dir = tempdir().unwrap();
    let gz = trim(dir.path(), &[]);
    assert_eq!(
        digest(&gz),
        SERIAL_GZ,
        "compressed bytes of the serial path changed; see this file's header \
         before re-pinning"
    );
}

#[test]
fn parallel_gzip_output_bytes_are_pinned() {
    let dir = tempdir().unwrap();
    let gz = trim(dir.path(), &["--cores", "2"]);
    assert_eq!(
        digest(&gz),
        PARALLEL_GZ,
        "compressed bytes of the --cores path changed; see this file's header \
         before re-pinning"
    );
}

#[test]
fn the_two_paths_differ_in_framing_but_not_in_payload() {
    // `--cores N` concatenates one gzip member per batch, so its framing differs
    // from the serial path by design. Pinned so the reason stays visible: equal
    // bytes here would mean the parallel path had stopped batching.
    assert_ne!(
        SERIAL_GZ, PARALLEL_GZ,
        "the parallel path emits its own member framing"
    );

    let serial = tempdir().unwrap();
    let parallel = tempdir().unwrap();
    let s = trim(serial.path(), &[]);
    let p = trim(parallel.path(), &["--cores", "2"]);

    // And the payload must not depend on --cores — the invariant CI checks through
    // `gzip -dc`, asserted here too so both halves live in one place.
    assert_eq!(
        gunzip(&s),
        gunzip(&p),
        "decompressed output must not depend on --cores"
    );
}

fn gunzip(path: &Path) -> Vec<u8> {
    use std::io::Read;
    let f = std::fs::File::open(path).expect("output exists");
    let mut d = flate2::read::MultiGzDecoder::new(f);
    let mut v = Vec::new();
    d.read_to_end(&mut v).expect("valid gzip");
    v
}
