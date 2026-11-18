//! Static assertions on `.github/workflows/release.yml`.
//!
//! The GA release-notes check has to sit in `check-release`, which publishes
//! nothing. `create-tag-and-release` pushes the tag before it creates the
//! release, so a check that failed there would leave a tag with no release
//! behind it — the v2.3.0 partial ship.
//!
//! Nothing can run after a failed job 0 as long as three things hold: every
//! other job needs `check-release`, no status function appears, and no
//! `continue-on-error` appears. All are properties of the file, so they are
//! asserted here rather than by a dispatch.
//!
//! A fourth property is not about job order: `Create release` reads the notes
//! path from a job-0 output, and an empty path is not refused by `gh` — it
//! posts an empty body and exits 0. So the link from the check to that flag is
//! pinned too.

use std::path::PathBuf;

/// The step name that validates the GA notes file.
const NOTES_STEP: &str = "Require curated GA release notes";

/// Every job the workflow is known to have.
///
/// Spelled out rather than counted: a new job has to be audited for `needs:`,
/// and failing here is how that audit gets requested.
const EXPECTED_JOBS: [&str; 9] = [
    "check-release",
    "build-binaries",
    "smoke-test-binaries",
    "docker-build",
    "docker-merge",
    "smoke-test-docker",
    "create-tag-and-release",
    "upload-binaries",
    "publish-crate",
];

/// Calls that would let a dependent job run after job 0 failed.
const STATUS_FUNCTIONS: [&str; 4] = ["always", "cancelled", "success", "failure"];

fn workflow() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(".github")
        .join("workflows")
        .join("release.yml");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

/// Whether `src` calls `name` as a function, with any spacing inside the call.
///
/// `always ()` is accepted by both GitHub and actionlint, so matching the
/// literal `always()` would miss it.
fn calls(src: &str, name: &str) -> bool {
    let mut from = 0;
    while let Some(i) = src[from..].find(name) {
        let start = from + i;
        let after = src[start + name.len()..].trim_start();
        if let Some(rest) = after.strip_prefix('(')
            && rest.trim_start().starts_with(')')
        {
            return true;
        }
        from = start + name.len();
    }
    false
}

/// Job name → its block, as `(name, lines)`, in file order.
///
/// Anchored to the `jobs:` mapping, because `on:` has two-space keys of its own
/// (`push:`) and treating one as a job would make every span wrong. A header is
/// any two-space key: GitHub allows `-`, `_` and digits in a job id, and
/// filtering on a narrower character set silently merges such a job into the
/// block above it, which makes the assertions below pass without reading it.
fn jobs(src: &str) -> Vec<(String, Vec<&str>)> {
    let lines: Vec<&str> = src.lines().collect();
    let start = lines
        .iter()
        .position(|l| *l == "jobs:")
        .expect("release.yml has no top-level `jobs:` key");

    let is_job_header = |l: &str| {
        l.starts_with("  ")
            && !l.starts_with("   ")
            && !l.trim_start().starts_with('#')
            && l.trim_end().ends_with(':')
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

/// A job's `needs:` value, inline (`[a, b]`) or block (`- a` lines).
fn needs_of(body: &[&str]) -> Option<String> {
    let at = body
        .iter()
        .position(|l| l.trim_start().starts_with("needs:"))?;
    let mut value = body[at]
        .trim_start()
        .trim_start_matches("needs:")
        .to_string();
    for line in &body[at + 1..] {
        let t = line.trim_start();
        if t.starts_with('-') {
            value.push(' ');
            value.push_str(t.trim_start_matches('-').trim());
        } else if !t.is_empty() && !t.starts_with('#') {
            break;
        }
    }
    Some(value)
}

#[test]
fn ga_notes_check_lives_in_check_release() {
    let src = workflow();
    assert!(
        src.contains(NOTES_STEP),
        "the GA release-notes check is gone from release.yml entirely"
    );

    for (name, body) in jobs(&src) {
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
fn release_workflow_cannot_run_past_a_failed_job_zero() {
    let src = workflow();
    for f in STATUS_FUNCTIONS {
        assert!(
            !calls(&src, f),
            "release.yml calls `{f}()`, which can let a job run after `check-release` failed"
        );
    }
    assert!(
        !src.contains("continue-on-error"),
        "release.yml uses `continue-on-error`: on a job it lets dependents run after a failure, \
         on a step it makes that step's `exit 1` inert"
    );
}

#[test]
fn every_other_job_needs_check_release() {
    let src = workflow();
    for (name, body) in jobs(&src) {
        if name == "check-release" {
            continue;
        }
        let needs = needs_of(&body).unwrap_or_else(|| {
            panic!("job `{name}` has no `needs:`, so a failed job 0 cannot stop it")
        });
        assert!(
            needs.contains("check-release"),
            "job `{name}` does not need `check-release`; its needs are:{needs}"
        );
    }
}

#[test]
fn the_notes_path_reaches_create_release() {
    // `gh release create --notes-file ""` is not refused — it posts an empty
    // body and exits 0 — and it would do so after the tag push. Each link in
    // the chain from the check to the flag is pinned.
    let src = workflow();
    for fragment in [
        "notes_path: ${{ steps.notes.outputs.notes_path }}",
        // Line-exact: `id: notes` is a prefix of `id: notes_check`.
        "\n        id: notes\n",
        "echo \"notes_path=$NOTES\" >> \"$GITHUB_OUTPUT\"",
        "--notes-file \"${{ needs.check-release.outputs.notes_path }}\"",
    ] {
        assert!(
            src.contains(fragment),
            "release.yml no longer contains `{fragment}`, so the GA body can reach \
             `--notes-file \"\"` after the tag push"
        );
    }

    // The check runs on exactly the dispatches that take the --notes-file arm.
    assert!(
        src.contains("steps.check.outputs.is_prerelease == 'false'"),
        "the notes check must be gated on the GA case, or a prerelease demands a notes file"
    );
    assert!(
        src.contains("--generate-notes --prerelease"),
        "the prerelease arm must keep --generate-notes"
    );
}

#[test]
fn job_parsing_sees_every_job() {
    // The assertions above are vacuous if this mis-parses: a job absorbed into
    // the block above it is never read, so it is never audited.
    let names: Vec<String> = jobs(&workflow()).into_iter().map(|(n, _)| n).collect();
    assert_eq!(
        names, EXPECTED_JOBS,
        "release.yml's jobs are not the set this test audits. If a job was added, add it to \
         EXPECTED_JOBS and make sure it needs `check-release`."
    );
}
