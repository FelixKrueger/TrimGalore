#!/bin/bash
# Real-data multi-member check. The bismark fixture I first reached for turned
# out to be a single member (`members=1`), so concatenated members were only
# covered by a synthetic 2-member file on the laptop. This builds a genuine
# 3-member archive out of real ENA reads and decodes it from a pipe.
set -u
BIN=/tmp/rgz_spike/crate/target/release/spike-nonseekable
ENA=/tmp/rgz_spike/SRR24827373_1.fastq.gz
cd /tmp/rgz_spike || exit 1
ms() { echo $(( ( $2 - $1 ) / 1000000 )); }

echo "=== exact decoded size of the ENA fixture (not gzip -l, which wraps at 4 GB) ==="
"$BIN" open "$ENA" > /dev/null 2> exact.txt
grep -oE 'OK bytes=[0-9]+ members=[0-9]+' exact.txt

if [ ! -s mm.gz ]; then
  echo "=== building a 3-member archive from real reads (3 GB of plain FASTQ) ==="
  gzip -dc "$ENA" | head -c 3000000000 > mm_plain
  split -n 3 -d mm_plain mm_part_
  for p in mm_part_00 mm_part_01 mm_part_02; do gzip -1 -c "$p" > "$p.gz"; done
  cat mm_part_00.gz mm_part_01.gz mm_part_02.gz > mm.gz
  ls -l mm_plain mm.gz | awk '{printf "  %.2f GB  %s\n", $5/1073741824, $9}'
fi

exp=$(md5sum mm_plain | cut -d' ' -f1)

echo
echo "=== decode the 3-member archive ==="
for mode in open pipe; do
  t0=$(date +%s%N)
  if [ "$mode" = open ]; then
    got=$(RGZ_THREADS=16 "$BIN" open mm.gz 2> mm_$mode.txt | md5sum | cut -d' ' -f1)
  else
    got=$(cat mm.gz | RGZ_THREADS=16 "$BIN" stream 2> mm_$mode.txt | md5sum | cut -d' ' -f1)
  fi
  t1=$(date +%s%N)
  tel=$(grep -oE 'mid_active=[0-9]+ mid_spawned=[0-9]+|members=[0-9]+' mm_$mode.txt | tr '\n' ' ')
  if [ "$got" = "$exp" ]; then
    echo "PASS  multimember_$mode  $(ms "$t0" "$t1") ms  [$tel]"
  else
    echo "FAIL  multimember_$mode  md5 mismatch: exp ${exp:0:12} got ${got:0:12}  $(head -1 mm_$mode.txt)"
  fi
done

echo
echo "=== control: md5 must reject a truncated payload ==="
if [ "$exp" = "$(head -c 1000 mm_plain | md5sum | cut -d' ' -f1)" ]; then
  echo "FAIL  md5_control"
else
  echo "PASS  md5_control — comparison can distinguish"
fi
