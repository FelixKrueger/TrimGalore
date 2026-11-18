# GA release notes

One file per GA release, named for its tag: `v2.6.0.md`. The release workflow reads it as the
release body. Prereleases are not covered — they keep `--generate-notes`.

## The mechanism, first

`check-release` — job 0 of the release workflow, which publishes nothing — fails a GA dispatch when
`.github/release-notes/<version>.md` is missing, empty, or only whitespace. The error names the path
to create.

That placement is the point. `create-tag-and-release` pushes the tag before it creates the release,
so a check that failed *there* would leave a tag on the remote with no release behind it — the shape
of the v2.3.0 partial ship. Validating in job 0 means that by the time anything is published the
file's existence is already proven. A CI test asserts the check stays in `check-release` and that the
workflow contains no status functions, either of which would let a later job run after job 0 failed.

Write the file on the release-prep branch, in the same commit as the version bump and the CHANGELOG
section. It is a trimmed restatement of that CHANGELOG section, so the material already exists.

## What a body contains

In order, nothing else.

**One line of framing.** What this release is. The CHANGELOG's own register — "The Clumpify Edition",
"The Formats Edition".

**Three to six highlights.** One sentence each, in the CHANGELOG's voice:

```
- **Headline** ([#NNN](https://github.com/FelixKrueger/TrimGalore/pull/NNN)) — one sentence.
```

The test for a highlight is whether a user would change what they do because of it. Dependency bumps,
CI, docs and internal refactors never qualify, however much work they were.

**A link to the CHANGELOG section** for the full record.

**Contributor acknowledgements, when the release has first-time contributors.** GitHub's own
heading and wording, so the credit reads the same as it always has:

```
## New Contributors
* @handle made their first contribution in https://github.com/FelixKrueger/TrimGalore/pull/NNN
```

This is the one component a curated body can lose by accident. `--generate-notes` appended it
automatically; `--notes-file` does not, so it has to be carried deliberately or the release page
simply stops crediting people.

Take the list from the auto-generated notes **at release time** — the preview GitHub shows when
drafting a release, or `gh api -X POST repos/FelixKrueger/TrimGalore/releases/generate-notes -f
tag_name=vX.Y.Z`. Do not regenerate it later to reconstruct an old release: the computation is not
stable over time. Regenerating v2.4.0's list in November 2026 returns two of the three contributors
its published body credits, because someone who was a first-time contributor then is not one now.
For past releases the published body is the only authoritative record.

**An install block.** Define it here rather than copying a previous release, because no previous
release body has one:

````
## Install

```bash
cargo install trim-galore
conda install -c bioconda trim-galore
docker pull ghcr.io/felixkrueger/trimgalore:latest
```
````

## What a body excludes

The full PR list, migration instructions, rationale, per-caller mappings. Migration instructions
belong in `docs/`; point at them in a clause rather than inlining them.

The one-sentence-per-entry rule that governs `CHANGELOG.md` governs this too. It is the same reader
asking the same question — *does this upgrade affect me?* — and a body that answers it in 16,000
characters has not answered it.

A finished body lands near 1,500–2,500 characters. That is a sanity range, not a gate.

## What these files are not

A record of what was announced, not a living document. Once a version ships, its file is not
rewritten. Editing the body on GitHub afterwards is fine and is not tracked back here — the file is
the source for *creation* only.
