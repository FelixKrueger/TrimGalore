#!/bin/bash
# Timing half of the oxy run. Separate script because `bc` is not installed on
# oxy, so the first driver produced empty durations; this uses integer
# nanosecond arithmetic in bash instead.
set -u
BIN=/tmp/rgz_spike/crate/target/release/spike-nonseekable
ENA=/tmp/rgz_spike/SRR24827373_1.fastq.gz
cd /tmp/rgz_spike || exit 1

ms() { echo $(( ( $2 - $1 ) / 1000000 )); }

echo "input: $(ls -l "$ENA" | awk '{printf "%.2f GB", $5/1073741824}') compressed"
echo "plain size: $(gzip -l "$ENA" 2>/dev/null | awk 'NR==2{print $2}') bytes (gzip -l, may wrap at 4 GB)"
echo

# Two runs each; report both so a cold-cache outlier is visible rather than averaged away.
timed() { # $1 label  $2 mode  $3 threads
  local label="$1" mode="$2" th="$3" t0 t1
  for rep in 1 2; do
    t0=$(date +%s%N)
    if [ "$mode" = open ]; then
      RGZ_THREADS="$th" "$BIN" open "$ENA" > /dev/null 2> "t_${label}.txt"
    else
      cat "$ENA" | RGZ_THREADS="$th" "$BIN" stream > /dev/null 2> "t_${label}.txt"
    fi
    t1=$(date +%s%N)
    printf "  %-22s rep%d  %6d ms   [%s]\n" "$label" "$rep" "$(ms "$t0" "$t1")" \
      "$(grep -oE 'mid_active=[0-9]+ mid_spawned=[0-9]+' "t_${label}.txt" | head -1)"
  done
}

echo "=== regular file via open() — parallel path ==="
timed file_t1  open 1
timed file_t4  open 4
timed file_t16 open 16
timed file_t32 open 32

echo
echo "=== same bytes down a pipe via stream_reader() — sequential by design ==="
timed pipe_t1  stream 1
timed pipe_t32 stream 32

echo
echo "=== references ==="
for rep in 1 2; do
  t0=$(date +%s%N); gzip -dc "$ENA" > /dev/null; t1=$(date +%s%N)
  printf "  %-22s rep%d  %6d ms\n" "gzip -dc" "$rep" "$(ms "$t0" "$t1")"
done
for rep in 1 2; do
  t0=$(date +%s%N); cat "$ENA" > /dev/null; t1=$(date +%s%N)
  printf "  %-22s rep%d  %6d ms\n" "cat (I/O floor)" "$rep" "$(ms "$t0" "$t1")"
done
