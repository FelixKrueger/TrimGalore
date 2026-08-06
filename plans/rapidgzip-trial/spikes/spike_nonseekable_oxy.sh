#!/bin/bash
# oxy driver for spike_nonseekable: real data, real core count, and a thread
# sweep so the parallel path is actually engaged — which the laptop run failed
# to do (30 copies of one text compressed to 1.7 MB and never left sequential).
set -u
BIN=/tmp/rgz_spike/crate/target/release/spike-nonseekable
cd /tmp/rgz_spike || exit 1

ENA=/tmp/rgz_spike/SRR24827373_1.fastq.gz
TGOUT=/home/fkrueger/bismark_benchmarks/RRBS_PE/synth_barcode/trim_out/synth_barcode_R1_val_1.fq.gz
BCL=/home/fkrueger/bismark_benchmarks/RRBS_PE/synth_bclconvert/synth_bclconvert_R1.fastq.gz

now() { date +%s.%N; }
el()  { echo "scale=2; ($2 - $1)/1" | bc; }

expected_md5() { # $1 path, $2 cache-key — gzip is the reference decoder
  local c="exp_$2.md5"
  if [ ! -s "$c" ]; then gzip -dc "$1" | md5sum | cut -d' ' -f1 > "$c"; fi
  cat "$c"
}

echo "=== fixtures ==="
for p in "$ENA" "$TGOUT" "$BCL"; do
  [ -s "$p" ] && ls -l "$p" | awk '{printf "  %.2f GB  %s\n", $5/1073741824, $9}' || echo "  MISSING $p"
done
echo "cores: $(nproc) logical"
echo

pass=0; fail=0
run_case() { # $1 label  $2 fixture  $3 key  $4 mode(open|pipe)  $5 threads
  local label="$1" fx="$2" key="$3" mode="$4" th="$5"
  local exp; exp=$(expected_md5 "$fx" "$key")
  local t0 t1 got meta="meta_${label}.txt"
  t0=$(now)
  if [ "$mode" = open ]; then
    got=$(RGZ_THREADS="$th" "$BIN" open "$fx" 2> "$meta" | md5sum | cut -d' ' -f1)
  else
    got=$(cat "$fx" | RGZ_THREADS="$th" "$BIN" stream 2> "$meta" | md5sum | cut -d' ' -f1)
  fi
  t1=$(now)
  local secs; secs=$(el "$t0" "$t1")
  local tel; tel=$(grep -oE 'mid_active=[0-9]+ mid_spawned=[0-9]+' "$meta" | head -1)
  local mem; mem=$(grep -oE 'members=[0-9]+' "$meta" | head -1)
  if [ "$got" = "$exp" ]; then
    echo "PASS  ${label}  ${secs}s  [$tel $mem]"; pass=$((pass+1))
  else
    echo "FAIL  ${label}  ${secs}s  md5 mismatch (exp ${exp:0:12} got ${got:0:12})  $(head -1 "$meta")"; fail=$((fail+1))
  fi
}

echo "=== correctness on real data: regular file vs pipe, 32 decoder threads ==="
run_case ena_file_t32          "$ENA"   ena  open  32
run_case ena_pipe_t32          "$ENA"   ena  pipe  32
run_case tgout_file_t32        "$TGOUT" tg   open  32
run_case tgout_pipe_t32        "$TGOUT" tg   pipe  32
run_case bclconvert_file_t32   "$BCL"   bcl  open  32
run_case bclconvert_pipe_t32   "$BCL"   bcl  pipe  32

echo
echo "=== does the parallel path engage on real data? thread sweep, ENA file, open() ==="
for th in 1 4 16 32 64; do run_case "ena_file_t${th}" "$ENA" ena open "$th"; done

echo
echo "=== the same sweep down a pipe — should not scale (sequential by design) ==="
for th in 1 32; do run_case "ena_pipe_t${th}" "$ENA" ena pipe "$th"; done

echo
echo "=== reference: gzip -dc wall time (single-threaded baseline) ==="
t0=$(now); gzip -dc "$ENA" > /dev/null; t1=$(now); echo "  gzip -dc: $(el "$t0" "$t1")s"

echo
echo "=== negative controls ==="
head -c 100000000 "$ENA" > /tmp/rgz_spike/trunc.gz
if cat /tmp/rgz_spike/trunc.gz | "$BIN" stream > /dev/null 2> meta_trunc.txt; then
  echo "FAIL  truncated_real_stream — accepted, exit 0"; fail=$((fail+1))
else
  echo "PASS  truncated_real_stream — $(grep -oE '(READ|OPEN|FINISH)_ERROR.*' meta_trunc.txt | head -1 | cut -c1-70)"; pass=$((pass+1))
fi
# Prove the md5 comparison can fail.
if [ "$(expected_md5 "$ENA" ena)" = "$(gzip -dc "$ENA" | head -c 1000 | md5sum | cut -d' ' -f1)" ]; then
  echo "FAIL  md5_control — a truncated payload compared equal"; fail=$((fail+1))
else
  echo "PASS  md5_control — md5 distinguishes full from truncated payload"; pass=$((pass+1))
fi

echo
echo "=== $pass passed, $fail failed ==="
