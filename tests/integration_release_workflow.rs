//! Static assertions on `.github/workflows/release.yml`.
//!
//! The GA release-notes check has to sit in `check-release`, which publishes
//! nothing. `create-tag-and-release` pushes the tag before it creates the
//! release, so a check that failed there would leave a tag with no release
//! behind it — the v2.3.0 partial ship.
//!
//! Nothing can run after a failed job 0 as long as two things hold: every other
//! job needs `check-release`, and the workflow has no status functions. Both are
//! properties of the file, so they are asserted here rather than by a dispatch.

use std::path::PathBuf;

/// The step name Step 2 of `plans/11112026_release-notes-policy/` installs.
const NOTES_STEP: &str = "Require curated GA release notes";

/// `always()` and friends would let a dependent job run after job 0 failed.
const STATUS_FUNCTIONS: [&str; 4] = ["always()", "cancelled()", "success()", "failure()"];

fn workflow() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(".github")
        .join("workflows")
        .join("release.yml");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Job name → its block, as `(name, lines)`, in file order.
///
/// Anchored to the `jobs:` mapping: `on:` has two-space keys of its own
/// (`push:`), and treating one as a job would make every span wrong.
fn jobs(src: &str) -> Vec<(String, Vec<&str>)> {
    let lines: Vec<&str> = src.lines().collect();
    let start = lines
        .iter()
        .position(|l| *l == "jobs:")
        .expect("release.yml has no top-level `jobs:` key");

    let is_job_header = |l: &str| {
        l.starts_with("  ")
            && !l.starts_with("   ")
            && l.trim_end().ends_with(':')
            && l.trim()
                .trim_end_matches(':')
                .chars()
                .all(|c| c.is_ascii_lowercase() || c == '-')
    };

    let headers: Vec<usize> = (start + 1..lines.len())
        .filter(|i| is_job_header(lines[*i]))
        .collect();

    headers
        .iter()
        .enumerate()
        .map(|(n, &i)| {
            let end = headers.get(n + 1).copied().unwrap_or(lines.len());
            let name = lines[i].trim().trim_end_matches(':').to_string();
            (name, lines[i..end].to_vec())
        })
        .collect()
}

#[test]
fn ga_notes_check_lives_in_check_release() {
    let src = workflow();
    assert!(
        src.contains(NOTES_STEP),
        "the GA release-notes check is gone from release.yml entirely"
    );

    let jobs = jobs(&src);
    assert!(jobs.len() > 1, "job parsing found {} jobs", jobs.len());

    for (name, body) in &jobs {
        let here = body.iter().any(|l| l.contains(NOTES_STEP));
        if name == "check-release" {
            assert!(
                here,
                "the GA release-notes check must be in `check-release`, which publishes nothing"
            );
        } else {
            assert!(
                !here,
                "the GA release-notes check is in `{name}`; a failure there can land after the \
                 tag push, which is the v2.3.0 partial ship"
            );
        }
    }
}

#[test]
fn release_workflow_has_no_status_functions() {
    let src = workflow();
    for f in STATUS_FUNCTIONS {
        assert!(
            !src.contains(f),
            "release.yml uses `{f}`, which can let a job run after `check-release` failed; \
             the no-status-function property is what makes the job-0 check sufficient"
        );
    }
}

#[test]
fn every_other_job_needs_check_release() {
    let src = workflow();
    for (name, body) in jobs(&src) {
        if name == "check-release" {
            continue;
        }
        let needs = body
            .iter()
            .find(|l| l.trim_start().starts_with("needs:"))
            .unwrap_or_else(|| {
                panic!("job `{name}` has no `needs:`, so a failed job 0 cannot stop it")
            });
        assert!(
            needs.contains("check-release"),
            "job `{name}` does not need `check-release`: {}",
            needs.trim()
        );
    }
}

#[test]
fn job_parsing_sees_the_real_jobs() {
    // The two assertions above are vacuous if this mis-parses, and `on:` has
    // two-space keys that look like job headers.
    let names: Vec<String> = jobs(&workflow()).into_iter().map(|(n, _)| n).collect();
    assert!(
        names.contains(&"check-release".to_string()),
        "parsed jobs {names:?} do not include check-release"
    );
    assert!(
        names.contains(&"create-tag-and-release".to_string()),
        "parsed jobs {names:?} do not include create-tag-and-release"
    );
    assert!(
        !names.contains(&"push".to_string()),
        "`on:` key parsed as a job: {names:?}"
    );
}
