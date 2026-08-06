#!/bin/bash
# Driver for spike_nonseekable. Builds gzip fixtures, runs the (format x source)
# matrix, and byte-compares decoded output against the original.
#
# Every invocation is bounded by `perl alarm` (alarm(2) survives execve), because
# a FIFO read that blocks would otherwise hang indefinitely and there is no
# timeout(1) on Darwin.
set -u
BIN="$1"      # path to the built spike binary
W="$2"        # scratch working directory
mkdir -p "$W" || exit 1
cd "$W" || exit 1

pass=0; fail=0
bounded() { perl -e 'alarm 60; exec @ARGV' "$@"; }

# ── fixtures ────────────────────────────────────────────────────────────────
# ~2.4 MB uncompressed: larger than the 64 KiB pipe buffer and any input page,
# so the incremental path is exercised rather than one convenient read.
if [ ! -s orig.fastq ]; then
  for i in $(seq 1 20000); do
    printf '@read%d\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$i"
  done > orig.fastq
fi
gzip -c orig.fastq > single.gz
head -c 500000 orig.fastq > part_a.txt; tail -c +500001 orig.fastq > part_b.txt
gzip -c part_a.txt > a.gz; gzip -c part_b.txt > b.gz
cat a.gz b.gz > multi.gz                 # concatenated members
bgzip -c orig.fastq > bgzf.gz            # BGZF framing + 28-byte EOF member
head -c 40000 single.gz > trunc.gz       # truncated: must error
printf 'sentinel\n' | cat orig.fastq - > orig_sentinel.fastq

echo "fixtures: orig=$(wc -c < orig.fastq)B single=$(wc -c < single.gz)B multi=$(wc -c < multi.gz)B bgzf=$(wc -c < bgzf.gz)B"
echo

# ── helpers ─────────────────────────────────────────────────────────────────
# $1 label  $2 expected-plain  $3.. command
check_identical() {
  local label="$1" expected="$2"; shift 2
  local out="out_${label}.bin" meta="meta_${label}.txt"
  "$@" > "$out" 2> "$meta"; local rc=$?
  if [ $rc -ne 0 ]; then
    echo "FAIL  $label — exit $rc: $(head -1 "$meta")"; fail=$((fail+1)); return
  fi
  if cmp -s "$out" "$expected"; then
    echo "PASS  $label — byte-identical ($(wc -c < "$out")B)  [$(grep -o 'members=[0-9]*' "$meta" || true) $(grep -o 'spawned_workers=[0-9]*' "$meta" || true)]"
    pass=$((pass+1))
  else
    echo "FAIL  $label — decoded output differs from the original"; fail=$((fail+1))
  fi
}

# $1 label  $2.. command   (expects a non-zero exit)
check_errors() {
  local label="$1"; shift
  local meta="meta_${label}.txt"
  "$@" > /dev/null 2> "$meta"; local rc=$?
  if [ $rc -ne 0 ]; then
    echo "PASS  $label — refused as expected: $(grep -oE '(READ|OPEN|FINISH)_ERROR.*' "$meta" | head -1 | cut -c1-90)"
    pass=$((pass+1))
  else
    echo "FAIL  $label — accepted truncated input, exit 0"; fail=$((fail+1))
  fi
}

# ── the matrix ──────────────────────────────────────────────────────────────
echo "── regular file (control: parallel path) ──"
check_identical "regular_single"  orig.fastq  bounded "$BIN" open single.gz
check_identical "regular_bgzf"    orig.fastq  bounded "$BIN" open bgzf.gz

echo
echo "── stdin from a pipe (stream_reader) ──"
sh -c "cat single.gz | perl -e 'alarm 60; exec @ARGV' '$BIN' stream" > out_pipe_single.bin 2> meta_pipe_single.txt
if cmp -s out_pipe_single.bin orig.fastq; then
  echo "PASS  pipe_single — byte-identical  [$(grep -o 'spawned_workers=[0-9]*' meta_pipe_single.txt)]"; pass=$((pass+1))
else echo "FAIL  pipe_single"; fail=$((fail+1)); fi

sh -c "cat multi.gz | perl -e 'alarm 60; exec @ARGV' '$BIN' stream" > out_pipe_multi.bin 2> meta_pipe_multi.txt
if cmp -s out_pipe_multi.bin orig.fastq; then
  echo "PASS  pipe_multi (concatenated members) — byte-identical  [$(grep -o 'members=[0-9]*' meta_pipe_multi.txt)]"; pass=$((pass+1))
else echo "FAIL  pipe_multi"; fail=$((fail+1)); fi

sh -c "cat bgzf.gz | perl -e 'alarm 60; exec @ARGV' '$BIN' stream" > out_pipe_bgzf.bin 2> meta_pipe_bgzf.txt
if cmp -s out_pipe_bgzf.bin orig.fastq; then
  echo "PASS  pipe_bgzf (BGZF + EOF member) — byte-identical  [$(grep -o 'members=[0-9]*' meta_pipe_bgzf.txt)]"; pass=$((pass+1))
else echo "FAIL  pipe_bgzf"; fail=$((fail+1)); fi

echo
echo "── named FIFO via open() — the auto-routing case ──"
F=fifo_$$
mkfifo "$F"
( cat single.gz > "$F" ) &
writer=$!
check_identical "fifo_open" orig.fastq bounded "$BIN" open "$F"
kill "$writer" 2>/dev/null; wait "$writer" 2>/dev/null

echo
echo "── the shapes TrimGalore is asked about ──"
sh -c "perl -e 'alarm 60; exec @ARGV' '$BIN' open /dev/stdin < single.gz" > out_devstdin.bin 2> meta_devstdin.txt
if cmp -s out_devstdin.bin orig.fastq; then
  echo "PASS  open_/dev/stdin_from_file — byte-identical  [$(grep -o 'spawned_workers=[0-9]*' meta_devstdin.txt)]"; pass=$((pass+1))
else echo "FAIL  open_/dev/stdin_from_file: $(head -1 meta_devstdin.txt)"; fail=$((fail+1)); fi

sh -c "cat single.gz | perl -e 'alarm 60; exec @ARGV' '$BIN' open /dev/stdin" > out_devstdin_pipe.bin 2> meta_devstdin_pipe.txt
if cmp -s out_devstdin_pipe.bin orig.fastq; then
  echo "PASS  open_/dev/stdin_from_pipe — byte-identical  [$(grep -o 'spawned_workers=[0-9]*' meta_devstdin_pipe.txt)]"; pass=$((pass+1))
else echo "FAIL  open_/dev/stdin_from_pipe: $(head -1 meta_devstdin_pipe.txt)"; fail=$((fail+1)); fi

/bin/bash -c "perl -e 'alarm 60; exec @ARGV' '$BIN' open <(cat single.gz)" > out_procsub.bin 2> meta_procsub.txt
if cmp -s out_procsub.bin orig.fastq; then
  echo "PASS  open_process_substitution — byte-identical  [$(grep -o 'spawned_workers=[0-9]*' meta_procsub.txt)]"; pass=$((pass+1))
else echo "FAIL  open_process_substitution: $(head -1 meta_procsub.txt)"; fail=$((fail+1)); fi

echo
echo "── negative controls ──"
sh -c "cat trunc.gz | perl -e 'alarm 60; exec @ARGV' '$BIN' stream" > /dev/null 2> meta_trunc.txt
if [ -s meta_trunc.txt ] && grep -qE '(READ|OPEN|FINISH)_ERROR' meta_trunc.txt; then
  echo "PASS  truncated_stream_errors — $(grep -oE '(READ|OPEN|FINISH)_ERROR.*' meta_trunc.txt | head -1 | cut -c1-80)"; pass=$((pass+1))
else
  echo "FAIL  truncated_stream_errors — no error reported"; fail=$((fail+1)); fi

# Prove cmp can fail, so every PASS above means something.
if cmp -s out_pipe_single.bin orig_sentinel.fastq; then
  echo "FAIL  cmp_control — a 9-byte difference was not detected"; fail=$((fail+1))
else
  echo "PASS  cmp_control — cmp detects a 9-byte difference"; pass=$((pass+1)); fi

echo
echo "=== $pass passed, $fail failed ==="
[ "$fail" -eq 0 ]
