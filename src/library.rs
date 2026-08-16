//! Library-prep presets (`--library <name>`), issue
//! [#440](https://github.com/FelixKrueger/TrimGalore/issues/440).
//!
//! A preset is nothing more than four end-clipping numbers under a kit name. The
//! membership criterion is deliberately narrow: a prep belongs here only if it
//! expands *purely* into 5'/3' clipping. That is why `--rrbs` stays its own flag
//! (it also changes adapter-trim behaviour and composes with `--non_directional`),
//! why `--clock` and `--implicon` stay out (they encode UMIs into read IDs), and
//! why NuGEN Ovation RRBS is absent (diversity trimming is a different algorithm,
//! and it must run *without* `--rrbs`).
//!
//! The numbers match nf-core/methylseq's presets of the same names, so a run
//! driven through that pipeline and a run driven straight through Trim Galore
//! trim identically.
//!
//! Presets may change between releases (see `CHANGELOG.md`). That is safe because
//! both trimming reports record the expanded values, so an old report stays
//! self-describing and a past run is reproducible from its own record.

/// The four end-clipping values a preset expands into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PresetClips {
    pub clip_r1: usize,
    pub clip_r2: usize,
    pub three_prime_clip_r1: usize,
    pub three_prime_clip_r2: usize,
}

/// A library-prep preset selected with `--library`.
///
/// One enum rather than one boolean per kit: a prep is a single choice, and an
/// enum makes "two kits at once" unrepresentable instead of needing a conflict
/// rule for every pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum LibraryPreset {
    /// NEBNext Enzymatic Methyl-seq (EM-seq).
    #[value(name = "emseq", alias = "em_seq", alias = "em-seq")]
    EmSeq,
    /// Accel-NGS Methyl-seq, sold as Swift and now as IDT xGen Methyl-seq.
    /// Aliased because the vendor renamed twice and users reach for whichever
    /// name they bought it under.
    #[value(name = "accel", alias = "swift", alias = "xgen")]
    Accel,
    /// Zymo-Seq / Pico Methyl-Seq.
    #[value(name = "zymo")]
    Zymo,
    /// Single-cell BS-seq (scBS-seq).
    #[value(name = "scbs", alias = "single_cell", alias = "single-cell")]
    ScBs,
    /// Post-Bisulfite Adapter Tagging.
    #[value(name = "pbat")]
    Pbat,
}

impl LibraryPreset {
    /// The name printed in reports and log lines, whichever alias was typed.
    pub fn canonical_name(self) -> &'static str {
        match self {
            LibraryPreset::EmSeq => "emseq",
            LibraryPreset::Accel => "accel",
            LibraryPreset::Zymo => "zymo",
            LibraryPreset::ScBs => "scbs",
            LibraryPreset::Pbat => "pbat",
        }
    }

    /// The clipping this preset expands into.
    pub fn clips(self) -> PresetClips {
        let (c1, c2, t1, t2) = match self {
            LibraryPreset::EmSeq => (10, 10, 10, 10),
            LibraryPreset::Accel => (10, 15, 10, 10),
            LibraryPreset::Zymo => (10, 10, 10, 10),
            LibraryPreset::ScBs => (6, 6, 6, 6),
            LibraryPreset::Pbat => (8, 8, 8, 8),
        };
        PresetClips {
            clip_r1: c1,
            clip_r2: c2,
            three_prime_clip_r1: t1,
            three_prime_clip_r2: t2,
        }
    }
}

/// One clip flag given on the command line over a preset that also sets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClipOverride {
    /// Flag name as the user would type it, e.g. `--clip_R1`.
    pub flag: &'static str,
    pub preset_value: usize,
    pub user_value: usize,
}

/// The four clip values a user typed, before any preset is folded in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct UserClips {
    pub clip_r1: Option<usize>,
    pub clip_r2: Option<usize>,
    pub three_prime_clip_r1: Option<usize>,
    pub three_prime_clip_r2: Option<usize>,
}

/// Clipping in force for a run, plus the provenance the reports need.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedClips {
    /// `None` when no `--library` was given; the other fields are then just the
    /// user's own values passed through.
    pub preset: Option<LibraryPreset>,
    pub clip_r1: Option<usize>,
    pub clip_r2: Option<usize>,
    pub three_prime_clip_r1: Option<usize>,
    pub three_prime_clip_r2: Option<usize>,
    /// Flags the user set that a preset would otherwise have supplied.
    pub overrides: Vec<ClipOverride>,
}

impl ResolvedClips {
    /// The clipping in force, rendered as the command line that would reproduce
    /// it. Printed in the log and in both trimming reports, so a run stays
    /// reproducible from its own record even if the preset changes later.
    pub fn flag_summary(&self) -> String {
        let mut parts = Vec::new();
        for (flag, value) in [
            ("--clip_R1", self.clip_r1),
            ("--clip_R2", self.clip_r2),
            ("--three_prime_clip_R1", self.three_prime_clip_r1),
            ("--three_prime_clip_R2", self.three_prime_clip_r2),
        ] {
            if let Some(v) = value {
                parts.push(format!("{flag} {v}"));
            }
        }
        parts.join(" ")
    }

    /// True when any of the four values is set, from either source.
    pub fn any_set(&self) -> bool {
        self.clip_r1.is_some()
            || self.clip_r2.is_some()
            || self.three_prime_clip_r1.is_some()
            || self.three_prime_clip_r2.is_some()
    }
}

/// Fold a preset (if any) into the user's clip flags.
///
/// Explicit wins, and every displaced value is recorded so the log and the
/// reports can name both numbers. An error instead of an override would push the
/// user back to writing all four values by hand, which is the thing presets exist
/// to remove; it also matches what `--rrbs` already does with `--clip_R2`.
pub fn resolve(preset: Option<LibraryPreset>, user: UserClips) -> ResolvedClips {
    let Some(preset) = preset else {
        return ResolvedClips {
            preset: None,
            clip_r1: user.clip_r1,
            clip_r2: user.clip_r2,
            three_prime_clip_r1: user.three_prime_clip_r1,
            three_prime_clip_r2: user.three_prime_clip_r2,
            overrides: Vec::new(),
        };
    };

    let clips = preset.clips();
    let mut overrides = Vec::new();
    let mut pick =
        |flag: &'static str, user_value: Option<usize>, preset_value: usize| match user_value {
            Some(v) => {
                overrides.push(ClipOverride {
                    flag,
                    preset_value,
                    user_value: v,
                });
                Some(v)
            }
            None => Some(preset_value),
        };

    ResolvedClips {
        preset: Some(preset),
        clip_r1: pick("--clip_R1", user.clip_r1, clips.clip_r1),
        clip_r2: pick("--clip_R2", user.clip_r2, clips.clip_r2),
        three_prime_clip_r1: pick(
            "--three_prime_clip_R1",
            user.three_prime_clip_r1,
            clips.three_prime_clip_r1,
        ),
        three_prime_clip_r2: pick(
            "--three_prime_clip_R2",
            user.three_prime_clip_r2,
            clips.three_prime_clip_r2,
        ),
        overrides,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table in issue #440, verbatim. These numbers are the feature.
    #[test]
    fn preset_table_matches_the_agreed_numbers() {
        let expect = |p: LibraryPreset, c1, c2, t1, t2| {
            assert_eq!(
                p.clips(),
                PresetClips {
                    clip_r1: c1,
                    clip_r2: c2,
                    three_prime_clip_r1: t1,
                    three_prime_clip_r2: t2,
                },
                "{} preset",
                p.canonical_name()
            );
        };
        expect(LibraryPreset::EmSeq, 10, 10, 10, 10);
        expect(LibraryPreset::Accel, 10, 15, 10, 10);
        expect(LibraryPreset::Zymo, 10, 10, 10, 10);
        expect(LibraryPreset::ScBs, 6, 6, 6, 6);
        expect(LibraryPreset::Pbat, 8, 8, 8, 8);
    }

    #[test]
    fn preset_supplies_all_four_values_when_user_gave_none() {
        let r = resolve(Some(LibraryPreset::Accel), UserClips::default());
        assert_eq!(r.clip_r1, Some(10));
        assert_eq!(r.clip_r2, Some(15));
        assert_eq!(r.three_prime_clip_r1, Some(10));
        assert_eq!(r.three_prime_clip_r2, Some(10));
        assert!(r.overrides.is_empty());
    }

    #[test]
    fn explicit_clip_flag_wins_over_the_preset() {
        let r = resolve(
            Some(LibraryPreset::EmSeq),
            UserClips {
                clip_r1: Some(12),
                ..UserClips::default()
            },
        );
        assert_eq!(r.clip_r1, Some(12));
        assert_eq!(
            r.clip_r2,
            Some(10),
            "untouched values still come from preset"
        );
    }

    #[test]
    fn an_override_records_both_values() {
        let r = resolve(
            Some(LibraryPreset::EmSeq),
            UserClips {
                clip_r1: Some(12),
                ..UserClips::default()
            },
        );
        assert_eq!(
            r.overrides,
            vec![ClipOverride {
                flag: "--clip_R1",
                preset_value: 10,
                user_value: 12,
            }]
        );
    }

    #[test]
    fn without_a_preset_user_values_pass_through_untouched() {
        let r = resolve(
            None,
            UserClips {
                clip_r2: Some(4),
                ..UserClips::default()
            },
        );
        assert_eq!(r.preset, None);
        assert_eq!(r.clip_r1, None);
        assert_eq!(r.clip_r2, Some(4));
        assert!(r.overrides.is_empty());
    }
}
