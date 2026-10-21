#!/bin/bash
# =============================================================================
# 05_deseq.sh - Build the count matrix and run the legacy DESeq test
# =============================================================================
# Identical to the original 2018 commands:
#   paste <(awk '{print $1"\t"$7"\t"$8"\t"$9}' WT_base_mRNA180209) \
#         <(awk '{print $7"\t"$8"\t"$9}' pol4ab_base_mRNA180209) > count table
#   Rscript DESeq.R WT_pol4ab_base.count        (bin/DESeq.R, unmodified)
# The DESeq step uses the correct replicate-wise counts in BOTH script sets
# (erroneous_script / corrected_script) - the historical FPKM bug never
# affected the differential-expression statistics.
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT/featurecount_DEG"

paste <(awk '{print $1"\t"$7"\t"$8"\t"$9}' WT_base_mRNA180209) \
      <(awk '{print $7"\t"$8"\t"$9}' pol4ab_base_mRNA180209) > WT_pol4ab_base.count

"$ROOT/envs/deseq/bin/Rscript" "$ROOT/bin/DESeq.R" WT_pol4ab_base.count

# restrict to protein-coding genes (39,962 IDs, list recovered from the
# original output WT_pol4_PCG.diff; the original ricev7.PCG.bed copies are
# dead symlinks - only the gene-ID column ever entered the results)
awk 'NR==FNR{p[$1]=1;next} p[$1]' "$ROOT/annotation/ricev7.PCG.bed" \
    WT_pol4.diff > WT_pol4_PCG.diff

echo "DESeq done. Significant PCGs (foldChange>=2 & padj<=0.01):"
awk '{if($5!="NA" && $5>=2 && $7<=0.01)print}' WT_pol4_PCG.diff | wc -l
echo "Significant PCGs (foldChange<=0.5 & padj<=0.01):"
awk '{if($5!="NA" && $5<=0.5 && $7<=0.01)print}' WT_pol4_PCG.diff | wc -l
