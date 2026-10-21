#!/bin/bash
# =============================================================================
# 04_featureCounts.sh - Fragment counting with featureCounts (subread 1.6.0)
# =============================================================================
# Identical parameters to the original 2018 command (bin/bsub_featureCounts.sh):
#   featureCounts -O -P -s 2 -p -d 10 -D 2000 -C -T 4 -a osativa7_newpol4ab.gtf \
#                  -o <name> rep1.bam rep2.bam rep3.bam
# One run per genotype (3 BAMs counted jointly), output named after the first
# 3 underscore fields of the BAM path, as in the original.
# The two header lines of the featureCounts output are stripped afterwards,
# matching the file format found in the original directory.
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/featurecount_DEG" && cd "$ROOT/featurecount_DEG"
export PATH="$ROOT/envs/subread/bin:$PATH"
cp "$ROOT/annotation/osativa7_newpol4ab.gtf" .

bash "$ROOT/bin/bsub_featureCounts.sh" \
  ../alignment/WT_base_mRNA180209_rep1/WT_base_mRNA180209_rep1.bam \
  ../alignment/WT_base_mRNA180209_rep2/WT_base_mRNA180209_rep2.bam \
  ../alignment/WT_base_mRNA180209_rep3/WT_base_mRNA180209_rep3.bam

bash "$ROOT/bin/bsub_featureCounts.sh" \
  ../alignment/pol4ab_base_mRNA180209_rep1/pol4ab_base_mRNA180209_rep1.bam \
  ../alignment/pol4ab_base_mRNA180209_rep2/pol4ab_base_mRNA180209_rep2.bam \
  ../alignment/pol4ab_base_mRNA180209_rep3/pol4ab_base_mRNA180209_rep3.bam

# strip the 2 header lines (original directory keeps headerless files)
for f in WT_base_mRNA180209 pol4ab_base_mRNA180209; do
  head -2 "$f" > "$f.header"; tail -n +3 "$f" > "$f.new" && mv "$f.new" "$f"
done
echo "featureCounts done."
