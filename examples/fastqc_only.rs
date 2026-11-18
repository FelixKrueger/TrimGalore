// Diagnostic-only: invoke the bundled fastqc-rust on a given path.
// Used to compare bundled vs Java FastQC output directly. Not shipped.
use std::path::Path;

fn main() -> anyhow::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: fastqc_only <input.fastq.gz> <outdir>");
        std::process::exit(1);
    }
    trim_galore::fastqc::run(
        Path::new(&args[1]),
        None,
        Some(Path::new(&args[2])),
        1,
        false,
        // Sniff the format from the extension, as the trimmed-output sites do.
        None,
    )?;
    Ok(())
}
