use std::env;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn git_short_hash() -> String {
    Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

fn build_epoch() -> u64 {
    match env::var("SOURCE_DATE_EPOCH") {
        // .trim() absorbs whitespace from shell expansions like
        // `SOURCE_DATE_EPOCH="$(git log -1 --format=%ct) "`; sign/format
        // strictness (no leading `+`, no underscores, no hex) is preserved.
        Ok(s) => s.trim().parse::<u64>().unwrap_or_else(|_| {
            panic!(
                "SOURCE_DATE_EPOCH must be a non-negative decimal seconds-since-epoch integer, got {s:?}"
            )
        }),
        Err(_) => SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before 1970-01-01")
            .as_secs(),
    }
}

fn format_iso8601_utc(epoch: u64) -> String {
    let secs_of_day = epoch % 86_400;
    let days = epoch / 86_400;
    let hour = secs_of_day / 3600;
    let minute = (secs_of_day % 3600) / 60;
    let second = secs_of_day % 60;

    let (year, month, day) = civil_from_days(days as i64);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

// Howard Hinnant's days-from-civil algorithm (public domain).
// Input: days since 1970-01-01.  Output: (year, month, day) on the proleptic Gregorian calendar.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    let year = y + i64::from(m <= 2);
    (year, m, d)
}

fn main() {
    let hash = git_short_hash();
    let timestamp = format_iso8601_utc(build_epoch());
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_else(|_| "unknown".to_string());
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_else(|_| "unknown".to_string());

    let version_body = format!("{hash} — {target_os}/{target_arch} — built {timestamp}");

    println!("cargo:rustc-env=GIT_SHORT_HASH={hash}");
    println!("cargo:rustc-env=BUILD_TIMESTAMP={timestamp}");
    println!("cargo:rustc-env=VERSION_BODY={version_body}");

    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    println!("cargo:rerun-if-env-changed=SOURCE_DATE_EPOCH");
}
