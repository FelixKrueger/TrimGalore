# Measurements — re-open behaviour and FIFO open semantics

Taken during manual review of `PLAN.md`, 2026-08-04, kernel Darwin 25.5.0 / macOS 26.5.1, arm64, rustc from `~/.cargo/bin`. **Every row below is Darwin**; the plan's Linux claims are reasoned, not measured, and PLAN §9 V9 makes CI measure the one that matters.
Probe sources: `probe.rs` (two opens + byte compare) and `probe2.rs` (fstat of the first handle), reproduced at the end.

Both probes mirror `format::detect_input_format`: open, `read` up to 4 bytes, then — **still holding the first handle**, as `format.rs:231` does — open again. `probe.rs` runs the second open on a thread with `recv_timeout(3s)` so a blocking open is reported rather than hanging.

## 1. Does re-opening restart? (`probe.rs`)

| # | input | open1 | open2 | restarts |
|---|---|---|---|---|
| A | regular file | `@MS2` | `@MS2` | **true** |
| B | `/dev/stdin` redirected from a regular file | `@MS2` | `_rea` | false |
| C | `/dev/stdin` from a pipe | `@MS2` | `_rea` | false |
| D | `<(cat reg.fq)` process substitution (bash) | `@MS2` | `_rea` | false |
| E | named FIFO, **writer exited** after writing | `@MS2` | *blocked >3 s* | **HANG** |
| F | named FIFO, writer still open | `@MS2` | `_rea` | false |

A–D reproduce PLAN §2.2 exactly, including the seekable-but-non-restarting case B. A is the positive control: the comparison can return `true`.

**E is new.** POSIX `open(O_RDONLY)` on a FIFO blocks until a writer opens; an existing reader holding the fd does not satisfy it. So the probe's own second open blocks indefinitely when the writer has exited — which is the common shape (`mkfifo f; cat small.fq > f &`).

E vs F is decided by whether the writer is still alive at second-open time, i.e. by payload size against the pipe buffer (64 KiB). A test using a small fixture is a scheduler coin-flip.

## 2. Can the file type be named without a second open? (`probe2.rs`)

| # | input | `File::metadata()` (fstat) | `fs::metadata(path)` |
|---|---|---|---|
| A | regular file | regular | regular |
| B | `/dev/stdin` from a regular file | **regular** | regular |
| C | `/dev/stdin` from a pipe | FIFO/pipe | FIFO/pipe |
| D | process substitution | FIFO/pipe | FIFO/pipe |
| E | named FIFO | FIFO/pipe | FIFO/pipe |

fstat on the handle already open in `detect_input_format` names every pipe and FIFO (C, D, E) — the three cases where a second open either consumes stream bytes or blocks. It reports **regular** for B, confirming PLAN §2.2: a type check alone cannot catch the shared-offset case, so it cannot replace the byte probe.

`std::os::unix::fs::FileTypeExt::is_fifo()` is std — no new dependency.

## 3. Consequences for the plan

1. Reject `is_fifo()` / `is_socket()` on the first handle **before** attempting a second open. Removes case E's hang, avoids consuming stream bytes, and leaves the byte probe to do the only job a type check cannot (case B).
2. A3/V5's risk is **flakiness**, not a deterministic hang: make the FIFO test's writer outlive the second open (or exceed the pipe buffer) and give the test an in-process timeout. `timeout(1)` is absent on Darwin and on GitHub's `macos-latest` leg, and the CI matrix is `[ubuntu-latest, macos-latest]`.

## 4. Probe sources

```rust
// probe.rs — two opens, byte compare, second open on a timed thread
use std::fs::File; use std::io::Read; use std::sync::mpsc; use std::thread;
use std::time::Duration;
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let mut f1 = File::open(&path).unwrap();
    let mut b1 = [0u8; 4];
    let n1 = f1.read(&mut b1).unwrap();
    let p2 = path.clone();
    let (tx, rx) = mpsc::channel();
    thread::spawn(move || {
        let r = (|| -> std::io::Result<(usize, [u8; 4])> {
            let mut g = File::open(&p2)?;
            let mut b = [0u8; 4];
            let n = g.read(&mut b)?;
            Ok((n, b))
        })();
        let _ = tx.send(r);
    });
    match rx.recv_timeout(Duration::from_secs(3)) {
        Ok(Ok((n2, b2))) => println!("restarts={}", b1[..n1] == b2[..n2]),
        Ok(Err(e)) => println!("second-open-failed: {e}"),
        Err(_) => println!("HANG"),
    }
    drop(f1);
    std::process::exit(0);
}
```

```rust
// probe2.rs — file type of the already-open handle vs the path
use std::fs::File; use std::os::unix::fs::FileTypeExt;
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let f = File::open(&path).unwrap();
    let ft = f.metadata().unwrap().file_type();
    println!("fifo={} sock={} chardev={} regular={}",
             ft.is_fifo(), ft.is_socket(), ft.is_char_device(), ft.is_file());
}
```

Harness notes: `nice(5) failed` lines in the FIFO runs are the sandbox refusing to renice a background job; they do not affect the result. Each FIFO used a fresh `$$`-suffixed name rather than being cleared.
