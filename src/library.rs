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
//!
//! `ClipFlag` also lives here. It names the four clip flags and the read each
//! belongs to, which the inert-flag warning and the `--rename` uBAM guard need
//! with no preset involved.

/// One of the four end-clipping flags.
///
/// An enum rather than a string table: the read a flag belongs to is one match
/// arm rather than a fact each consumer re-derives, and a fifth flag will not
/// compile until every consumer handles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClipFlag {
    ClipR1,
    ClipR2,
    ThreePrimeClipR1,
    ThreePrimeClipR2,
}

impl ClipFlag {
    /// The closed set, for `any_effective` and the drift tests. This order is
    /// deliberately neither display order, so no display site can iterate it and
    /// still look correct.
    pub const ALL: [ClipFlag; 4] = [
        ClipFlag::ClipR2,
        ClipFlag::ClipR1,
        ClipFlag::ThreePrimeClipR2,
        ClipFlag::ThreePrimeClipR1,
    ];

    /// Flag name as the user would type it. These four strings reach the JSON
    /// report, so they are output.
    pub fn name(self) -> &'static str {
        match self {
            ClipFlag::ClipR1 => "--clip_R1",
            ClipFlag::ClipR2 => "--clip_R2",
            ClipFlag::ThreePrimeClipR1 => "--three_prime_clip_R1",
            ClipFlag::ThreePrimeClipR2 => "--three_prime_clip_R2",
        }
    }

    /// A Read 2 flag.
    pub fn read_2(self) -> bool {
        matches!(self, ClipFlag::ClipR2 | ClipFlag::ThreePrimeClipR2)
    }

    /// This flag names a read the run processes. `paired` means `--paired`; the
    /// modes that pair two files without it read no clip flag at all.
    pub fn applies(self, paired: bool) -> bool {
        paired || !self.read_2()
    }
}

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
    pub flag: ClipFlag,
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
    /// The clipping in force, rendered as the command line that would reproduce it.
    /// R2 values are listed for paired runs only.
    pub fn flag_summary(&self, paired: bool) -> String {
        let mut parts = Vec::new();
        // This order is user-visible in both trimming reports: 5' before 3',
        // Read 1 before Read 2.
        for flag in [
            ClipFlag::ClipR1,
            ClipFlag::ClipR2,
            ClipFlag::ThreePrimeClipR1,
            ClipFlag::ThreePrimeClipR2,
        ] {
            if !flag.applies(paired) {
                continue;
            }
            if let Some(v) = self.value(flag) {
                parts.push(format!("{} {v}", flag.name()));
            }
        }
        parts.join(" ")
    }

    /// The value in force for `flag`, preset folded in.
    pub fn value(&self, flag: ClipFlag) -> Option<usize> {
        match flag {
            ClipFlag::ClipR1 => self.clip_r1,
            ClipFlag::ClipR2 => self.clip_r2,
            ClipFlag::ThreePrimeClipR1 => self.three_prime_clip_r1,
            ClipFlag::ThreePrimeClipR2 => self.three_prime_clip_r2,
        }
    }

    /// True when a clip value in force belongs to a read this run processes.
    /// Says nothing about whether the *mode* honours clip flags; the caller's
    /// hardtrim terms cover that.
    pub fn any_effective(&self, paired: bool) -> bool {
        ClipFlag::ALL
            .iter()
            .any(|&f| f.applies(paired) && self.value(f).is_some())
    }
}

impl UserClips {
    /// The value the user typed for `flag`, before any preset.
    pub fn value(&self, flag: ClipFlag) -> Option<usize> {
        match flag {
            ClipFlag::ClipR1 => self.clip_r1,
            ClipFlag::ClipR2 => self.clip_r2,
            ClipFlag::ThreePrimeClipR1 => self.three_prime_clip_r1,
            ClipFlag::ThreePrimeClipR2 => self.three_prime_clip_r2,
        }
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
    let mut pick = |flag: ClipFlag, user_value: Option<usize>, preset_value: usize| match user_value
    {
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
        clip_r1: pick(ClipFlag::ClipR1, user.clip_r1, clips.clip_r1),
        clip_r2: pick(ClipFlag::ClipR2, user.clip_r2, clips.clip_r2),
        three_prime_clip_r1: pick(
            ClipFlag::ThreePrimeClipR1,
            user.three_prime_clip_r1,
            clips.three_prime_clip_r1,
        ),
        three_prime_clip_r2: pick(
            ClipFlag::ThreePrimeClipR2,
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
                flag: ClipFlag::ClipR1,
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

    /// `read_2` is a match on the variant, not a parse of the name, so the name
    /// serves here as an independent oracle.
    #[test]
    fn every_clip_flag_knows_which_read_it_belongs_to() {
        for f in ClipFlag::ALL {
            assert_eq!(
                f.read_2(),
                f.name().ends_with("_R2"),
                "{} has read_2 = {}",
                f.name(),
                f.read_2()
            );
        }
    }

    /// A fifth variant makes the match non-exhaustive, which forces a visit here.
    /// Nothing forces it into `ALL`, whose length is a literal.
    #[test]
    fn all_holds_every_clip_flag_once() {
        for f in ClipFlag::ALL {
            match f {
                ClipFlag::ClipR1
                | ClipFlag::ClipR2
                | ClipFlag::ThreePrimeClipR1
                | ClipFlag::ThreePrimeClipR2 => {}
            }
        }
        let mut names: Vec<&str> = ClipFlag::ALL.iter().map(|f| f.name()).collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), 4, "ALL must hold each flag exactly once");
    }

    /// These four strings reach the JSON report, so a typo is a silent output
    /// change. The `_R2` oracle above passes with or without the leading dashes.
    #[test]
    fn clip_flag_names_are_the_command_line_spellings() {
        assert_eq!(ClipFlag::ClipR1.name(), "--clip_R1");
        assert_eq!(ClipFlag::ClipR2.name(), "--clip_R2");
        assert_eq!(ClipFlag::ThreePrimeClipR1.name(), "--three_prime_clip_R1");
        assert_eq!(ClipFlag::ThreePrimeClipR2.name(), "--three_prime_clip_R2");
    }

    /// #450 — a Read 2 value reaches no clip site on a single-end run.
    #[test]
    fn any_effective_discounts_read_2_on_single_end() {
        let r2_only = resolve(
            None,
            UserClips {
                clip_r2: Some(5),
                ..UserClips::default()
            },
        );
        assert!(!r2_only.any_effective(false));
        assert!(r2_only.any_effective(true));

        let tp_r2_only = resolve(
            None,
            UserClips {
                three_prime_clip_r2: Some(7),
                ..UserClips::default()
            },
        );
        assert!(!tp_r2_only.any_effective(false));
        assert!(tp_r2_only.any_effective(true));

        // Read 1 counts on either shape, and a live R1 flag is not masked by an
        // inert R2 one.
        for user in [
            UserClips {
                clip_r1: Some(3),
                ..UserClips::default()
            },
            UserClips {
                three_prime_clip_r1: Some(3),
                ..UserClips::default()
            },
            UserClips {
                clip_r1: Some(3),
                clip_r2: Some(5),
                ..UserClips::default()
            },
        ] {
            let r = resolve(None, user);
            assert!(r.any_effective(false));
            assert!(r.any_effective(true));
        }

        assert!(!resolve(None, UserClips::default()).any_effective(true));
    }

    #[test]
    fn flag_summary_lists_all_four_for_paired() {
        let r = resolve(Some(LibraryPreset::Accel), UserClips::default());
        assert_eq!(
            r.flag_summary(true),
            "--clip_R1 10 --clip_R2 15 --three_prime_clip_R1 10 --three_prime_clip_R2 10"
        );
    }

    #[test]
    fn flag_summary_omits_read_2_for_single_end() {
        let r = resolve(Some(LibraryPreset::Accel), UserClips::default());
        assert_eq!(
            r.flag_summary(false),
            "--clip_R1 10 --three_prime_clip_R1 10"
        );
    }
}
