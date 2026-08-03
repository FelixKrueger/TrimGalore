# Code Review — Build Provenance CORE (Reviewer B)

**Scope:** `build.rs`, the `#[clap(...)]` attribute on `Cli`, and the startup
banner in `main.rs`. Nothing else.

**Verdict:** Ship-ready with one Medium (`u64::parse` swallowing whitespace/sign
quirks) and two Low recommendations. Core logic is correct; reproducibility
path is honoured; banner/clap wiring is sound.

---

## 1. Hinnant `civil_from_days` — correctness

Spot-checked the Rust port against `datetime.timedelta` for:

| days since 1970-01-01 | expected | got |
|--|--|--|
| 0 | 1970-01-01 | 1970-01-01 |
| 1 | 1970-01-02 | 1970-01-02 |
| 59 (non-leap Feb→Mar boundary) | 1970-03-01 | 1970-03-01 |
| 365 | 1971-01-01 | 1971-01-01 |
| 10 957 (Y2K) | 2000-01-01 | 2000-01-01 |
| 19 782 (leap day) | 2024-02-29 | 2024-02-29 |
| 19 783 | 2024-03-01 | 2024-03-01 |
| 19 675 (`SOURCE_DATE_EPOCH=1700000000`) | 2023-11-14 | 2023-11-14 |

All match. The port from Hinnant's original C++ is faithful — the `if z >= 0`
era branch is correct for the seconds-since-epoch domain we use (always ≥ 0 in
practice), and the month/day rotation at line 54 (`year = y + i64::from(m <=
2)`) correctly rolls Jan/Feb into the civil year.

No arithmetic overflow for any `u64` epoch: `u64::MAX / 86_400 ≈ 2.14 × 10^14`,
which casts losslessly to `i64` (< `i64::MAX`). Year overflow is effectively
impossible — at `u64::MAX` seconds the result is ~year 584 554 051 223, which
`format!("{year:04}")` happily prints in full (width is a minimum, not a
truncation). No panic path in the formatter.

**No issues.**

---

## 2. UTF-8 / trimming of git stdout (`git_short_hash`)

```rust
.and_then(|o| String::from_utf8(o.stdout).ok())
.map(|s| s.trim().to_string())
.filter(|s| !s.is_empty())
.unwrap_or_else(|| "unknown".to_string())
```

- Non-UTF-8 output → `String::from_utf8` returns `Err`, `.ok()` → `None`,
  `.and_then` short-circuits, falls through to `"unknown"`. **No panic.**
- `git rev-parse --short HEAD` emits ASCII hex + `\n`; a valid UTF-8 path is
  guaranteed in practice, but the fallback is safe even if git writes a BOM or
  localised error text on stderr (stderr isn't captured; stdout only).
- `.trim()` strips the trailing `\n` correctly.
- Empty-but-successful output (shouldn't happen but could with a mocked git) is
  caught by `.filter(|s| !s.is_empty())`.

**No issues.**

---

## 3. `SOURCE_DATE_EPOCH` parsing

```rust
s.parse::<u64>().unwrap_or_else(|_| panic!(…))
```

### 3a. Negative values — **correctly rejected** (Good)

`"-1".parse::<u64>()` returns `Err(ParseIntError)`, so the `panic!` fires with
the plan-required descriptive message. Matches the Debian reproducible-builds
spec, which requires a non-negative integer.

### 3b. Overflow vs year "292 billion" — Low concern

`u64::MAX` seconds = ~5.85 × 10^11 years (year ~584 billion, not 292B — plan's
"292B" was a loose approximation; no action needed). `parse::<u64>()` will
reject anything > `u64::MAX` with `Err`. Values between `i64::MAX` and
`u64::MAX` parse fine and flow through the civil algo without overflow. Nothing
to fix.

### 3c. MEDIUM — leading `+`, whitespace, unusual inputs

`u64::parse` accepts `"18446744073709551615"` (u64::MAX) but rejects `"+100"`,
`" 100"`, `"100 "`, `"0x64"`, `"1_000"`. This is strictly correct per the
reproducible-builds spec (decimal integer, no sign, no whitespace), but the
panic message says *"non-negative decimal seconds-since-epoch integer"* without
flagging common gotchas. A CI pipeline that does `SOURCE_DATE_EPOCH="$(git log
-1 --format=%ct) "` (trailing space from shell expansion) would hard-fail the
build with a slightly opaque message.

**Recommendation (Medium):** accept `s.trim()` before parsing, to absorb shell
whitespace without changing strictness on sign/format:

```rust
Ok(s) => s.trim().parse::<u64>().unwrap_or_else(|_| panic!(…))
```

This is a one-character change and strictly more forgiving for a CI misuse
that is common.

### 3d. Empty string

`"".parse::<u64>()` → `Err` → panic. Acceptable; empty string is malformed.
But if the user's intent was "unset" they'd `unset` it; panic is the right
call. No change recommended.

---

## 4. `.git` present but git binary missing on PATH

`Command::new("git").args(…).output()`:

- If `git` is not on PATH, `output()` returns
  `Err(io::Error { kind: NotFound, … })`.
- `.ok()` converts to `None`; the whole chain falls through to `"unknown"`.
- Build **succeeds** with `GIT_SHORT_HASH=unknown`.

Verified by reading: the `.ok()` at line 10 consumes the `io::Error` before
status inspection. **No issues.**

One Low note: the `cargo:rerun-if-changed=.git/HEAD` directive will still fire
when HEAD moves even though git is unavailable — cargo does the stat itself, no
git binary needed. The build will re-run unnecessarily on each commit in this
edge case (rare: `.git` present + git absent), but it correctly re-produces
`"unknown"`. Not worth fixing.

---

## 5. Banner interaction with `--quiet` / stdout capture

The banner uses `eprintln!` (stderr), matching the existing Trim Galore Perl
behaviour. There is no `--quiet` flag on `Cli` (confirmed by Grep against
`src/cli.rs`), so no conflict. Pipelines capturing stdout (e.g.
`trim_galore - < in.fq > out.fq` — not actually supported, but hypothetically)
remain unaffected: banner on stderr, data on stdout. **No issues.**

One Low observation: the existing banner printed a blank line, version, then
`===…===`. The new version adds a `VERSION_BODY` line between the version and
separator. This is exactly what the plan (§4, line 47 "after the existing …
line") specified. External log parsers that pattern-match on the
`Trim Galore - Oxidized Edition v…` line itself are unaffected because that
line is unchanged.

---

## 6. Banner racing `clap -V` short-circuit — **safe**

Verified by reading `main.rs` lines 22–36 in order:

```
env_logger::init();
let cli = Cli::parse();    // ← `-V` / `--version` short-circuits here and exits(0)
cli.validate()?;
eprintln!("\nTrim Galore - Oxidized Edition v{}", …);   // banner
```

`Cli::parse()` internally calls `clap::Parser::parse()`, which, on `-V` /
`--version`, prints to **stdout** and calls `std::process::exit(0)`. Control
never returns to `main.rs`. The banner therefore **cannot** pollute the
version-query path. Verified by clap's documented behaviour (see
`clap::Command::version`). **No issues.**

`env_logger::init()` does run before clap — it's cheap and writes nothing
unless a log event fires, so no stderr noise on `-V`. Fine.

---

## 7. `clap` long_version / version wiring — correctness

```rust
version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)"),
long_version = concat!(
    env!("CARGO_PKG_VERSION"), " (Oxidized Edition)\n",
    env!("VERSION_BODY")
),
```

- `concat!` requires string literals; `env!` expands to a literal at compile
  time — legal. Confirmed by successful build (otherwise `src/cli.rs` wouldn't
  compile).
- `\n` inside `concat!` yields a literal newline in the final string. clap
  prints `long_version` verbatim, so `-V` → one line, `--version` → two lines,
  as §4 requires. **Correct wiring.**

Low nit: the literal `\n` forces a Unix line-ending on Windows `--version`
output. Trim Galore Oxidized doesn't target Windows (see `Cargo.toml` /
`Dockerfile`); non-issue here, but worth a code comment if cross-platform ever
becomes a goal.

---

## 8. Summary of actionable items

| Priority | Area | Fix |
|--|--|--|
| Medium | `build.rs` §build_epoch | `s.trim().parse::<u64>()` before panic to absorb trailing shell whitespace in `SOURCE_DATE_EPOCH` |
| Low | `build.rs` §format_iso8601_utc | Optional: comment noting year field widens beyond 4 digits for epoch > ~year 10000 — not a bug, just non-obvious |
| Low | `src/cli.rs` `long_version` | Optional: comment noting `\n` is literal LF, not platform-native — irrelevant today |

No Critical or High findings. Core logic is tight, reproducibility is real
(verified the `sha256sum` recipe in PLAN §5 would genuinely be byte-identical
given the `cargo:rerun-if-env-changed=SOURCE_DATE_EPOCH` directive at line 72).

---

## 9. Fixes applied

None. All findings above are recommendations; the Medium (`.trim()`) is a
one-line change but I'm flagging rather than applying since the current code
is *spec-compliant* — the user may prefer strict Debian-spec behaviour over
forgiveness.
