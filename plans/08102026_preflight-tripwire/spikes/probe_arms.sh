#!/bin/zsh
# Probe: what does each dispatch arm actually create?
# Snapshot dir recursively before/after; print the created set.
set -u
TG="$1"; ROOT="$2"; REPO="$3"
rm -rf "$ROOT"; mkdir -p "$ROOT"

mk_fastq() {  # $1=path $2=prefix
  mkdir -p "$(dirname "$1")"
  : > "$1"
  for i in {0..39}; do
    print -r -- "@$2_read$i" >> "$1"
    print -r -- "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" >> "$1"
    print -r -- "+" >> "$1"
    print -r -- "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> "$1"
  done
}

snap() { find "$1" -type f | sed "s|^$1/||" | sort }

arm() {  # $1=name  rest=args ; fixtures already staged by caller in $D
  local name="$1"; shift
  local D="$ROOT/$name"
  local before after
  before=$(snap "$D")
  ( cd "$D" && "$TG" "$@" > stdout.log 2> stderr.log )
  local rc=$?
  after=$(snap "$D")
  print -r -- "########## $name (exit=$rc) :: $*"
  comm -13 <(print -r -- "$before") <(print -r -- "$after") | grep -v '^std\(out\|err\)\.log$'
  if [[ $rc -ne 0 ]]; then print -r -- "--- stderr tail ---"; tail -4 "$D/stderr.log"; fi
  print -r -- ""
}

# ---- 1. SE FASTQ trim, plain
D=$ROOT/se_trim; mk_fastq "$D/in/a.fastq" A; arm se_trim in/a.fastq

# ---- 2. SE FASTQ trim + fastqc
D=$ROOT/se_trim_fastqc; mk_fastq "$D/in/a.fastq" A; arm se_trim_fastqc --fastqc in/a.fastq

# ---- 3. SE FASTQ trim, two inputs, -o
D=$ROOT/se_trim_o; mk_fastq "$D/in/a.fastq" A; mk_fastq "$D/in/b.fastq" B
mkdir -p "$ROOT/se_trim_o/out"; arm se_trim_o -o out in/a.fastq in/b.fastq

# ---- 4. paired FASTQ trim, R1/R2 in DIFFERENT dirs (#398 rule)
D=$ROOT/pe_trim; mk_fastq "$D/r1/s_R1.fastq" R1; mk_fastq "$D/r2/s_R2.fastq" R2
arm pe_trim --paired r1/s_R1.fastq r2/s_R2.fastq

# ---- 5. paired + retain_unpaired
D=$ROOT/pe_retain; mk_fastq "$D/s_R1.fastq" R1; mk_fastq "$D/s_R2.fastq" R2
arm pe_retain --paired --retain_unpaired s_R1.fastq s_R2.fastq

# ---- 6. paired + passthrough
D=$ROOT/pe_pt; mk_fastq "$D/s_R1.fastq" R1; mk_fastq "$D/s_R2.fastq" R2; mk_fastq "$D/s_I1.fastq" I1
arm pe_pt --paired --passthrough s_I1.fastq s_R1.fastq s_R2.fastq

# ---- 7. paired + fastqc
D=$ROOT/pe_fastqc; mk_fastq "$D/s_R1.fastq" R1; mk_fastq "$D/s_R2.fastq" R2
arm pe_fastqc --paired --fastqc s_R1.fastq s_R2.fastq

# ---- 8. SE + demux
D=$ROOT/se_demux; mk_fastq "$D/a.fastq" A
print -r -- "sampleA\tACGT" > "$D/bc.txt"; print -r -- "sampleB\tTTTT" >> "$D/bc.txt"
arm se_demux --demux bc.txt a.fastq

# ---- 9. hardtrim5
D=$ROOT/ht5; mk_fastq "$D/in/a.fastq" A; arm ht5 --hardtrim5 20 in/a.fastq

# ---- 10. clock
D=$ROOT/clock; mk_fastq "$D/c_R1.fastq" R1; mk_fastq "$D/c_R2.fastq" R2
arm clock --clock --paired c_R1.fastq c_R2.fastq

# ---- 11. implicon
D=$ROOT/implicon; mk_fastq "$D/i_R1.fastq" R1; mk_fastq "$D/i_R2.fastq" R2
arm implicon --implicon --paired i_R1.fastq i_R2.fastq

# ---- 12. clump_only SE
D=$ROOT/co_se; mk_fastq "$D/a.fastq" A; arm co_se --clump_only a.fastq

# ---- 13. clump_only paired (R1/R2 different dirs)
D=$ROOT/co_pe; mk_fastq "$D/r1/s_R1.fastq" R1; mk_fastq "$D/r2/s_R2.fastq" R2
arm co_pe --clump_only --paired r1/s_R1.fastq r2/s_R2.fastq

# ---- 14. trim + clumpify (temp files?)
D=$ROOT/clumpify; mk_fastq "$D/a.fastq" A; arm clumpify --clumpify a.fastq

# ---- 15. SE uBAM out (FASTQ in)
D=$ROOT/ubam_se; mk_fastq "$D/a.fastq" A; arm ubam_se --output-format ubam a.fastq

# ---- 16. paired uBAM out, two FASTQ in
D=$ROOT/ubam_pe2; mk_fastq "$D/r1/s_R1.fastq" R1; mk_fastq "$D/r2/s_R2.fastq" R2
arm ubam_pe2 --paired --output-format ubam r1/s_R1.fastq r2/s_R2.fastq

# ---- 17. paired uBAM out from single interleaved BAM (ORPHAN: no pre-flight)
D=$ROOT/ubam_pe1; mkdir -p "$D"; cat "$REPO/test_files/ubam_paired_test.bam" > "$D/ip.bam"
arm ubam_pe1 --paired --output-format ubam ip.bam

# ---- 18. paired FASTQ out from single interleaved BAM (ORPHAN: no pre-flight)
D=$ROOT/fq_pe1; mkdir -p "$D"; cat "$REPO/test_files/ubam_paired_test.bam" > "$D/ip.bam"
arm fq_pe1 --paired ip.bam

# ---- 19. SE trim from uBAM input, FASTQ out
D=$ROOT/bamin_se; mkdir -p "$D"; cat "$REPO/test_files/ubam_test.bam" > "$D/s.bam"
arm bamin_se s.bam

# ---- 20. clump_only uBAM SE
D=$ROOT/co_ubam_se; mkdir -p "$D"; cat "$REPO/test_files/ubam_test.bam" > "$D/s.bam"
arm co_ubam_se --clump_only --output-format ubam s.bam

# ---- 21. clump_only uBAM paired shape B (single interleaved)
D=$ROOT/co_ubam_pe1; mkdir -p "$D"; cat "$REPO/test_files/ubam_paired_test.bam" > "$D/ip.bam"
arm co_ubam_pe1 --clump_only --paired --output-format ubam ip.bam

# ---- 22. clump_only uBAM paired shape A (two files)
D=$ROOT/co_ubam_pe2; mkdir -p "$D"; mk_fastq "$D/s_R1.fastq" R1; mk_fastq "$D/s_R2.fastq" R2
arm co_ubam_pe2 --clump_only --paired --output-format ubam s_R1.fastq s_R2.fastq

# ---- 23. hardtrim5 uBAM out
D=$ROOT/ht5_ubam; mk_fastq "$D/a.fastq" A; arm ht5_ubam --hardtrim5 20 --output-format ubam a.fastq

# ---- 24. no_report_file, SE
D=$ROOT/se_norep; mk_fastq "$D/a.fastq" A; arm se_norep --no_report_file a.fastq
