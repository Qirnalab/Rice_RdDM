#!/bin/bash
# =============================================================================
# 06_fpkm.sh  [CORRECTED VERSION]
# =============================================================================
# Same pipeline as erroneous_script/ steps 00-05 (identical alignment,
# counting and DESeq), with the FPKM formula CORRECTED:
#
#   erroneous (original 2018):  rep_i = 1000 * counts_rep1   / (library_size_i * length)
#   corrected (this script):    rep_i = 1000 * counts_rep_i  / (library_size_i * length)
#
# Each replicate's FPKM now uses that replicate's OWN fragment count
# (columns $7, $8, $9 of the featureCounts output = rep1, rep2, rep3).
#
# Library size = total fragments seen by featureCounts for that sample
# (sum over all status categories of the .summary file, in millions) -
# same denominator rule as the original.
#
# Worked example (gene LOC_Os01g01010, length 3,030 bp, original 2018 counts
# WT = 1112/1246/951):
#   erroneous  Rep1/2/3 = 12.1538/10.1939/14.6796, mean 12.3424  (published)
#   corrected  Rep1/2/3 = 12.1538/11.4223/12.5542, mean 12.0434  (revised table)
# Rep1 is identical in both (its numerator was always correct); only
# Rep2/Rep3 change, and they become genuine per-replicate measurements.
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
# CORRECTED FORMULA: rep_i uses counts column $7 (rep1), $8 (rep2), $9 (rep3)
# ---------------------------------------------------------------------------
awk -v L1="$W1" -v L2="$W2" -v L3="$W3" \
  '{rep1=1000*$7/(L1*$6);rep2=1000*$8/(L2*$6);rep3=1000*$9/(L3*$6);mean=(rep1+rep2+rep3)/3;print $1"\t"mean}' \
  WT_base_mRNA180209 > WT_base_mRNA180209.FPKM
awk -v L1="$P1" -v L2="$P2" -v L3="$P3" \
  '{rep1=1000*$7/(L1*$6);rep2=1000*$8/(L2*$6);rep3=1000*$9/(L3*$6);mean=(rep1+rep2+rep3)/3;print $1"\t"mean}' \
  pol4ab_base_mRNA180209 > pol4ab_base_mRNA180209.FPKM

# merge the two mean-FPKM columns onto the DESeq PCG table
paste WT_base_mRNA180209.FPKM <(cut -f2 pol4ab_base_mRNA180209.FPKM) \
  | perl "$ROOT/bin/map1.pl" - WT_pol4_PCG.diff > WT_pol4_PCG.diff_FPKM

# ---------------------------------------------------------------------------
# Replicate-level table in the format of the revised Supplemental Table 4
# (Gene, WT Rep1/2/3/Mean, osnrpd1-1 Rep1/2/3/Mean) - CORRECTED values,
# each replicate computed from its own count.
# ---------------------------------------------------------------------------
awk -v W1="$W1" -v W2="$W2" -v W3="$W3" -v P1="$P1" -v P2="$P2" -v P3="$P3" '
  NR==FNR {c1[$1]=$7; c2[$1]=$8; c3[$1]=$9; l[$1]=$6; next}
           {d1[$1]=$7; d2[$1]=$8; d3[$1]=$9}
  END {for (g in c1) {L=l[g];
     w1=1000*c1[g]/(W1*L);w2=1000*c2[g]/(W2*L);w3=1000*c3[g]/(W3*L);
     p1=1000*d1[g]/(P1*L);p2=1000*d2[g]/(P2*L);p3=1000*d3[g]/(P3*L);
     printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", \
            g,w1,w2,w3,(w1+w2+w3)/3,p1,p2,p3,(p1+p2+p3)/3}}' \
  WT_base_mRNA180209 pol4ab_base_mRNA180209 | sort -k1,1 \
  > FPKM_replicates_corrected.tsv
echo -e "Gene\tWT_Rep1\tWT_Rep2\tWT_Rep3\tWT_Mean\tosnrpd1-1_Rep1\tosnrpd1-1_Rep2\tosnrpd1-1_Rep3\tosnrpd1-1_Mean" \
  | cat - FPKM_replicates_corrected.tsv > tmp && mv tmp FPKM_replicates_corrected.tsv

echo "done. Outputs: WT_base_mRNA180209.FPKM, pol4ab_base_mRNA180209.FPKM, WT_pol4_PCG.diff_FPKM,"
echo "              FPKM_replicates_corrected.tsv (replicate-level table, revised-table format)"
