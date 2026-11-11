# Reconciliation measurement — threading.md memory table

Host: dockyard-oxy-0 (Intel Xeon 6975P-C, 128 cores, 991 GiB RAM), Linux/x86_64
Fixture: SRR24827373 Buckberry PE pair, 84M pairs, **65 bp** (verified 200k/200k reads both mates)
Invocation: `trim_galore --paired --cores N -o OUT R1.fastq.gz R2.fastq.gz` (gz in, gz out; matches scripts/benchmark.sh)
Method: /usr/bin/time -v peak RSS, 3 reps per cell, fixture staged on local NVMe
Guards: GNU-time capability probe; disk guard with an impossible-threshold positive control; per-run exit status, worker-count banner, and gz-output assertions all recorded

Published for comparison (threading.md:53-58, v2.1.0-beta.7): c4 = 73.0 MB, c8 = 91.6 MB

| build | ref | allocator | cores | max MiB | min MiB | mean MiB | spread MiB |
|---|---|---|---|--:|--:|--:|--:|
| beta7 | v2.1.0-beta.7 | stock | c4 | 84.3 | 74.1 | 78.8 | 10.2 |
| beta7 | v2.1.0-beta.7 | stock | c8 | 103.4 | 72.5 | 87.7 | 30.9 |
| v230 | c6528e5 (v2.3.0) | stock | c4 | 64.5 | 56.0 | 60.7 | 8.5 |
| v230 | c6528e5 (v2.3.0) | stock | c8 | 74.2 | 73.7 | 74.0 | 0.5 |
| dev | a107833 (dev) | mimalloc | c4 | 106.6 | 101.7 | 104.3 | 5.0 |
| dev | a107833 (dev) | mimalloc | c8 | 163.7 | 150.7 | 157.0 | 13.0 |

## Findings

1. **beta.7 reproduces the published table.** Published c4 = 73.0 MB and c8 = 91.6 MB both fall inside
   the measured beta.7 ranges (c4 74.1–84.3 MiB, c8 72.5–103.4 MiB), and wall time matches to ~2%
   (1:22.08 vs 1:20.75 at c4). The table is sound for the fixture it was measured on.
2. **The fixture is 65 bp, not 150 bp.** `threading.md:63` calls it "an 84M-read 150-bp fixture" and
   `:61` says "on these 150 bp reads". Both are wrong, and the mislabel is what makes the figure look
   applicable to 150 bp short-read data.
3. **v2.3.0 did not regress; current `dev` did.** On the identical fixture and host, peak RSS at
   `--cores 8` is 74.2 MiB on v2.3.0 and 163.7 MiB on `dev` — a **2.2x increase**, and 1.65x at c4.
   v2.3.0 is actually *lower* and far more stable than beta.7.
4. **Leading suspect: the mimalloc swap (#437, `2ab0918`), unconfirmed.** `dev` differs from v2.3.0 by
   the allocator *and* every commit since the tag, so the allocator is not yet isolated. #437's commit
   message reports a -24% peak-RSS move on the `--clumpify` path; the plain paired path moves the other
   way. mimalloc also buys wall time here (c8: 0:52 dev vs 1:07 v2.3.0), so this looks like a
   speed-for-memory trade rather than a defect.
   **Isolating experiment:** build `2ab0918` and `2ab0918^` and measure both at c8. One variable.
5. **Run-to-run spread is large enough to matter.** beta.7 varied 30.9 MiB across three identical c8
   runs. The published table is a single run per condition, so its 3-significant-figure precision
   overstates its reliability. v2.3.0 was the tightest (0.5 MiB).

## Status

Unreleased: v2.3.0 is what users install and what the nf-core container pins, so finding 3 affects
nobody today. It would ship with the next release.

## Isolation run — is mimalloc the cause?

Same host, same staged 65 bp fixture, `--paired`, gz in/out, 3 reps per cell.
Binary-level control: `strings | grep -ci mimalloc` = 7 for mimalloc builds, 0 for the others.

| build | ref | allocator | cores | mean MiB | max MiB | mean wall |
|---|---|---|---|--:|--:|--:|
| mimparent | `2ab0918^` | stock | c8 | 74.7 | 80.9 | 1:07.1 |
| mimchild | `2ab0918` | mimalloc | c8 | 164.3 | 167.7 | 0:51.8 |
| devnomi | `origin/dev + strip` | stock | c8 | 76.9 | 79.3 | 1:07.1 |
| dev | `origin/dev` | mimalloc | c8 | 146.0 | 148.6 | 0:52.8 |
| devnomi | `origin/dev + strip` | stock | c4 | 76.8 | 81.9 | 1:48.3 |
| dev | `origin/dev` | mimalloc | c4 | 97.6 | 104.7 | 1:42.1 |

### Verdict: mimalloc is the cause. Two independent one-variable tests agree.

- **`2ab0918^` -> `2ab0918`** (the allocator commit alone): 74.7 -> 164.3 MiB at c8 = **2.20x**.
- **`dev` minus `#[global_allocator]` -> `dev`** (the allocator alone, at today's tip):
  76.9 -> 146.0 MiB at c8 = **1.90x**; at c4 **1.27x**.
- **No code component.** `devnomi` at c8 (76.9 MiB) lands on v2.3.0's 74.0 MiB, so nothing
  in the ~40 commits since the tag raised plain-path memory. The allocator accounts for all of it.
- **It is a deliberate-looking trade, not a leak.** mimalloc buys ~28-30% wall at c8
  (1:07 -> 0:52) and ~6% at c4, for ~1.9-2.2x peak RSS. Page retention, not growth.

### Why it matters beyond the docs

1. Unreleased: v2.3.0 ships the stock allocator, so no user is affected yet.
2. #437's commit message reports mimalloc *reducing* `--clumpify` peak RSS by 24%. If both hold, the
   same swap moves the two paths in opposite directions, and the clumpify memory model's
   `STATIC_TRIM_BYTES = 224 MiB` was calibrated with mimalloc in place.
3. It compounds the #1903 failure mode: at 65 bp the plain path goes ~74 -> ~146 MiB, and 150 bp data
   is ~2.3x that again before any cgroup page-cache charge. A 1 GB budget stops looking generous.


## The c4 gap: a 2x2, and an allocator/reservation interaction

Bisected to `1595d59` (`perf(clump): cut --clumpify peak RSS 13-17%`), which replaced a flat
`reads.len() * 300` per-batch output-buffer reservation with a content-derived
`estimated_record_bytes` sum. That function is on the shared plain path, not just the clumpy one.

Mean peak RSS (MiB), 65 bp Buckberry pair, 3 reps per cell, same host:

| allocator | reservation | c4 | c8 | build |
|---|---|--:|--:|---|
| stock | flat | 59.2 | 74.0 | v2.3.0 as released |
| stock | estimator | 76.8 | 76.9 | dev, allocator stripped |
| mimalloc | flat | 104.0 | 161.0 | dev + reservation reverted |
| mimalloc | estimator | 97.6 | 146.0 | dev as it stands |

### The reservation change flips sign with the allocator

- Under the **stock** allocator it costs memory: +17.5 MiB at c4 and
  +3.0 MiB at c8. A constant 1,228,800-byte request every batch
  recycles perfectly on a free list; a per-batch-varying size does not.
- Under **mimalloc** it saves memory: -6.4 MiB at c4 and
  -15.0 MiB at c8. Size-class segregation handles the variance,
  so the padded flat request just retains ~40% more per buffer.

**So there is nothing to fix.** In the configuration that ships (mimalloc) the change is a modest win;
reverting it would cost ~15 MiB at c8. The c4 gap that prompted this only exists in the
stock-allocator-plus-estimator combination, which was never released and never will be -- it was an
artifact of holding the allocator constant to attribute the effect.

### Attribution of the v2.3.0 -> dev rise

At c8, mimalloc accounts for 96% of the move and at c4 54%, so the
threading docs' "mostly `mimalloc`" is accurate at both core counts. Because the two factors interact,
the split depends on the order of attribution; "most" is the robust claim, a precise percentage is not.

### Caveat

All of this is 65 bp reads on one host. The reservation effect could differ at 150 bp, where per-batch
variance is larger in absolute terms. Not measured.

