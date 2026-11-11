#!/bin/zsh
set -u
TG="$1"; ROOT="$2"; REPO="$3"
rm -rf "$ROOT"; mkdir -p "$ROOT"

mk_fastq() {  # $1=path $2=idprefix  -- SAME ids across mates when prefix identical
  mkdir -p "$(dirname "$1")"
  : > "$1"
  for i in {0..39}; do
    print -r -- "@${2}read$i" >> "$1"
    print -r -- "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" >> "$1"
    print -r -- "+" >> "$1"
    print -r -- "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> "$1"
  done
}
snap() { find "$1" \( -type f -o -type d \) | sed "s|^$1/*||" | sed '/^$/d' | sort }
arm() {
  local name="$1"; shift
  local D="$ROOT/$name"; local before after
  before=$(snap "$D")
  ( cd "$D" && "$TG" "$@" > stdout.log 2> stderr.log ); local rc=$?
  after=$(snap "$D")
  print -r -- "########## $name (exit=$rc) :: $*"
  comm -13 <(print -r -- "$before") <(print -r -- "$after") | grep -v '^std\(out\|err\)\.log$'
  [[ $rc -ne 0 ]] && { print -r -- "--- stderr ---"; tail -3 "$D/stderr.log" }
  print -r -- ""
}

# demux, with REAL tabs
D=$ROOT/demux; mk_fastq "$D/a.fastq" A
printf 'sampleA\tACGTACGT\nsampleB\tTTTTGGGG\n' > "$D/bc.txt"
arm demux --demux bc.txt a.fastq

# demux + fastqc
D=$ROOT/demux_fqc; mk_fastq "$D/a.fastq" A
printf 'sampleA\tACGTACGT\nsampleB\tTTTTGGGG\n' > "$D/bc.txt"
arm demux_fqc --demux bc.txt --fastqc a.fastq

# passthrough, ids in lockstep
D=$ROOT/pt; mk_fastq "$D/s_R1.fastq" S_; mk_fastq "$D/s_R2.fastq" S_; mk_fastq "$D/s_I1.fastq" S_
arm pt --paired --passthrough s_I1.fastq s_R1.fastq s_R2.fastq

# passthrough + fastqc
D=$ROOT/pt_fqc; mk_fastq "$D/s_R1.fastq" S_; mk_fastq "$D/s_R2.fastq" S_; mk_fastq "$D/s_I1.fastq" S_
arm pt_fqc --paired --fastqc --passthrough s_I1.fastq s_R1.fastq s_R2.fastq

# clumpify with cores (temp bins?)
D=$ROOT/clumpify; mk_fastq "$D/a.fastq" A; arm clumpify --clumpify --cores 2 a.fastq

# clumpify paired
D=$ROOT/clumpify_pe; mk_fastq "$D/s_R1.fastq" S_; mk_fastq "$D/s_R2.fastq" S_
arm clumpify_pe --paired --clumpify --cores 2 s_R1.fastq s_R2.fastq

# clump_only + fastqc
D=$ROOT/co_fqc; mk_fastq "$D/a.fastq" A; arm co_fqc --clump_only --fastqc a.fastq

# fastqc_args -o DIR : the documented residual
D=$ROOT/fqc_odir; mk_fastq "$D/a.fastq" A; mkdir -p "$D/qc"
arm fqc_odir --fastqc --fastqc_args "-o qc" a.fastq

# basename, paired
D=$ROOT/basename_pe; mk_fastq "$D/s_R1.fastq" S_; mk_fastq "$D/s_R2.fastq" S_
arm basename_pe --paired --basename fixed s_R1.fastq s_R2.fastq

# basename, SE two inputs (collapses -> should be refused)
D=$ROOT/basename_se; mk_fastq "$D/a.fastq" A; mk_fastq "$D/b.fastq" B
arm basename_se --basename fixed a.fastq b.fastq

# gzipped input -> gz outputs
D=$ROOT/gz; mk_fastq "$D/a.fastq" A; gzip -c "$D/a.fastq" > "$D/a.fastq.gz"; rm "$D/a.fastq"
arm gz a.fastq.gz

# multi-pair paired (2 pairs)
D=$ROOT/multipair; for p in 1 2; do mk_fastq "$D/p${p}_R1.fastq" S_; mk_fastq "$D/p${p}_R2.fastq" S_; done
arm multipair --paired p1_R1.fastq p1_R2.fastq p2_R1.fastq p2_R2.fastq

# hardtrim3 multi-input
D=$ROOT/ht3; mk_fastq "$D/a.fastq" A; mk_fastq "$D/b.fastq" B; arm ht3 --hardtrim3 20 a.fastq b.fastq

# clock multi-pair
D=$ROOT/clock2; for p in 1 2; do mk_fastq "$D/c${p}_R1.fastq" S_; mk_fastq "$D/c${p}_R2.fastq" S_; done
arm clock2 --clock --paired c1_R1.fastq c1_R2.fastq c2_R1.fastq c2_R2.fastq

# retain_unpaired that actually strands reads (short R2)
D=$ROOT/retain_real; mk_fastq "$D/s_R1.fastq" S_
: > "$D/s_R2.fastq"
for i in {0..39}; do
  print -r -- "@S_read$i" >> "$D/s_R2.fastq"; print -r -- "ACGTA" >> "$D/s_R2.fastq"
  print -r -- "+" >> "$D/s_R2.fastq"; print -r -- "IIIII" >> "$D/s_R2.fastq"
done
arm retain_real --paired --retain_unpaired s_R1.fastq s_R2.fastq

# uBAM out + fastqc (what extension does fastqc-rust see?)
D=$ROOT/ubam_fqc; mk_fastq "$D/a.fastq" A; arm ubam_fqc --output-format ubam --fastqc a.fastq

# uBAM out SE + no_report_file
D=$ROOT/ubam_norep; mk_fastq "$D/a.fastq" A; arm ubam_norep --output-format ubam --no_report_file a.fastq
