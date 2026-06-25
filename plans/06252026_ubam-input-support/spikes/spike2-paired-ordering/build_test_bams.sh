#!/usr/bin/env bash
# Build small synthetic uBAMs in three orderings to test samtools sort -n vs collate.
# Output: input.bam (interleaved), sorted_n.bam, collate.bam in same dir.
set -euo pipefail
cd "$(dirname "$0")"

# Generate a SAM with R1/R2 records for 5 read pairs, interleaved R1,R2,R1,R2,...
# Flag 77  = 0x4D = paired + unmapped + mate-unmapped + first-in-pair  (0x01|0x04|0x08|0x40)
# Flag 141 = 0x8D = paired + unmapped + mate-unmapped + second-in-pair (0x01|0x04|0x08|0x80)

{
  printf '@HD\tVN:1.6\tSO:unsorted\n'
  for i in 01 02 03 04 05 06 07 08 09 10; do
    printf 'read_%s\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' "$i"
    printf 'read_%s\t141\t*\t0\t0\t*\t*\t0\t0\tTGCATGCATG\tIIIIIIIIII\n' "$i"
  done
} > input.sam

samtools view -bS input.sam > input.bam

echo "=== input.bam (interleaved by construction): name, flag ==="
samtools view input.bam | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "=== samtools sort -n input.bam (queryname sort) ==="
samtools sort -n -o sorted_n.bam input.bam
samtools view sorted_n.bam | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "=== samtools collate input.bam (re-interleave for paired processing) ==="
samtools collate -O input.bam tmp_collate_prefix | samtools view - | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "=== Now construct an adversarial grouped input (all R1, then all R2), sort -n it ==="
{
  printf '@HD\tVN:1.6\tSO:unsorted\n'
  for i in 01 02 03 04 05 06 07 08 09 10; do
    printf 'read_%s\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n' "$i"
  done
  for i in 01 02 03 04 05 06 07 08 09 10; do
    printf 'read_%s\t141\t*\t0\t0\t*\t*\t0\t0\tTGCATGCATG\tIIIIIIIIII\n' "$i"
  done
} | samtools view -bS - > grouped.bam

echo "--- grouped.bam (R1 x10 then R2 x10): ---"
samtools view grouped.bam | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "--- samtools sort -n grouped.bam: ---"
samtools sort -n -o grouped_sorted_n.bam grouped.bam
samtools view grouped_sorted_n.bam | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "--- samtools collate grouped.bam: ---"
samtools collate -O grouped.bam tmp_collate_prefix2 | samtools view - | awk '{print NR": "$1, "flag="$2}'

echo ""
echo "=== Adversarial-2: queryname-tied R1s then queryname-tied R2s but with names that sort lexically same ==="
echo "(skipping — the names above already test the worst case: identical name R1 and R2 differ only by flag)"
