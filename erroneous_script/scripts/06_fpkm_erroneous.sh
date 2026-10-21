#!/bin/bash
# =============================================================================
# 06_fpkm.sh  [ERRONEOUS VERSION - reproduces the original published values]
# =============================================================================
# This script reproduces, verbatim, the FPKM calculation of the original 2018
# analysis INCLUDING ITS ERROR. It is provided so that reviewers can verify
# that the published replicate FPKM values follow from this formula.
#
# THE ERROR: the numerator of ALL THREE replicates is the read count of
# replicate 1 (column $7 of the featureCounts output). Only the library-size
# denominators differ between replicates. The correct formula (see
# corrected_script/) must use each replicate's own count ($7, $8, $9).
#
#   erroneous:  rep_i = 1000 * counts_rep1 / (library_size_i * length)
#   correct:    rep_i = 1000 * counts_rep_i / (library_size_i * length)
#
# Consequence: the three "replicate" FPKM values are not independent
# measurements - their ratios are fixed by the library sizes alone - and the
# reported mean is biased toward replicate 1 whenever replicate 1 deviates
# from the other replicates.
#
# Library size = total fragments seen by featureCounts for that sample
# (sum over all status categories of the .summary file, in millions).
# The 2018 run used hardcoded numbers; here they are computed from the
# current summaries by the same rule. Values obtained in our reproduction:
#   WT_base      : 30.288408 / 36.078303 / 25.022860   (2018: 30.196068/36.001589/25.000466)
#   pol4ab_base  : 28.519484 / 26.749251 / 25.721068   (2018: 28.513593/26.664223/25.622575)
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT/featurecount_DEG"

libsize () {  # $1 = sample prefix -> prints "L1 L2 L3" (millions)
  awk -v pre="$1" 'BEGIN{f=pre".summary"} NR>1{for(i=2;i<=NF;i++) c[i-1]+=$i}
                   END{printf "%.6f %.6f %.6f", c[1]/1e6, c[2]/1e6, c[3]/1e6}' "$1.summary"
}
read -r W1 W2 W3 <<< "$(libsize WT_base_mRNA180209)"
read -r P1 P2 P3 <<< "$(libsize pol4ab_base_mRNA180209)"
echo "library sizes (millions): WT $W1 $W2 $W3 | pol4ab $P1 $P2 $P3"

# ---------------------------------------------------------------------------
# THE ERRONEOUS FORMULA (identical structure to the original run.sh):
# NOTE $7 is used as the numerator in all three terms.
# ---------------------------------------------------------------------------
awk -v L1="$W1" -v L2="$W2" -v L3="$W3" \
  '{rep1=1000*$7/(L1*$6);rep2=1000*$7/(L2*$6);rep3=1000*$7/(L3*$6);mean=(rep1+rep2+rep3)/3;print $1"\t"mean}' \
  WT_base_mRNA180209 > WT_base_mRNA180209.FPKM
awk -v L1="$P1" -v L2="$P2" -v L3="$P3" \
  '{rep1=1000*$7/(L1*$6);rep2=1000*$7/(L2*$6);rep3=1000*$7/(L3*$6);mean=(rep1+rep2+rep3)/3;print $1"\t"mean}' \
  pol4ab_base_mRNA180209 > pol4ab_base_mRNA180209.FPKM

# merge the two mean-FPKM columns onto the DESeq PCG table
# (original command: paste ... | perl map1.pl - WT_pol4_PCG.diff)
paste WT_base_mRNA180209.FPKM <(cut -f2 pol4ab_base_mRNA180209.FPKM) \
  | perl "$ROOT/bin/map1.pl" - WT_pol4_PCG.diff > WT_pol4_PCG.diff_FPKM

# ---------------------------------------------------------------------------
# Replicate-level table in the format of published Supplemental Table 4
# (Gene, WT Rep1/2/3/Mean, osnrpd1-1 Rep1/2/3/Mean) - ERRONEOUS values.
# All three WT columns derive from the WT replicate-1 count only (same for
# osnrpd1-1): this reproduces the replicate columns of the published table.
# ---------------------------------------------------------------------------
awk -v W1="$W1" -v W2="$W2" -v W3="$W3" -v P1="$P1" -v P2="$P2" -v P3="$P3" '
  NR==FNR {c1[$1]=$7; l[$1]=$6; next}
           {c2[$1]=$7}
  END {for (g in c1) {L=l[g];
     w1=1000*c1[g]/(W1*L);w2=1000*c1[g]/(W2*L);w3=1000*c1[g]/(W3*L);
     p1=1000*c2[g]/(P1*L);p2=1000*c2[g]/(P2*L);p3=1000*c2[g]/(P3*L);
     printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", \
            g,w1,w2,w3,(w1+w2+w3)/3,p1,p2,p3,(p1+p2+p3)/3}}' \
  WT_base_mRNA180209 pol4ab_base_mRNA180209 | sort -k1,1 \
  > FPKM_replicates_erroneous.tsv
echo -e "Gene\tWT_Rep1\tWT_Rep2\tWT_Rep3\tWT_Mean\tosnrpd1-1_Rep1\tosnrpd1-1_Rep2\tosnrpd1-1_Rep3\tosnrpd1-1_Mean" \
  | cat - FPKM_replicates_erroneous.tsv > tmp && mv tmp FPKM_replicates_erroneous.tsv

echo "done. Outputs: WT_base_mRNA180209.FPKM, pol4ab_base_mRNA180209.FPKM, WT_pol4_PCG.diff_FPKM,"
echo "              FPKM_replicates_erroneous.tsv (replicate-level table, published-table format)"
