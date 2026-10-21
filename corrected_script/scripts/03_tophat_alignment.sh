#!/bin/bash
# =============================================================================
# 03_tophat_alignment.sh - TopHat2 alignment of the 6 libraries
# =============================================================================
# The tophat command (bin/rna-seq.sh) is byte-identical to the original 2018
# script except for the index path (original /Share path is gone):
#   tophat --bowtie1 -o $pre -G osativa7_newpol4ab.gtf --no-coverage-search \
#          -p 16 --library-type fr-firststrand <index> PE1 PE2
# Two run modes: SLURM (default; 6 jobs in parallel, as in the original
# cluster submission) or LOCAL (sequential).
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT/alignment"
export PATH="$ROOT/envs/tophat/bin:$PATH"
SAMPLES="WT_base_mRNA180209_rep1 WT_base_mRNA180209_rep2 WT_base_mRNA180209_rep3 \
pol4ab_base_mRNA180209_rep1 pol4ab_base_mRNA180209_rep2 pol4ab_base_mRNA180209_rep3"

MODE=${1:-slurm}
for S in $SAMPLES; do
  rm -rf "$S"
  if [ "$MODE" = slurm ] && command -v sbatch >/dev/null; then
    sbatch -J "${S:0:12}" -c 16 -o "../logs/${S}.slurm.out" -e "../logs/${S}.slurm.err" \
      --wrap="export PATH=$ROOT/envs/tophat/bin:\$PATH; bash rna-seq.sh ../rawdata/${S}_PE1.fq.gz ../rawdata/${S}_PE2.fq.gz"
  else
    bash rna-seq.sh "../rawdata/${S}_PE1.fq.gz" "../rawdata/${S}_PE2.fq.gz" > "../logs/${S}.tophat.log" 2>&1
  fi
done

# wait until all 6 samples produced their final BAM (~1.5 h in parallel)
while [ "$(ls */accepted_hits.bam 2>/dev/null | wc -l)" -lt 6 ]; do sleep 60; done
# also give the last jobs a minute to finish writing unmapped.bam / logs
sleep 120

# original convention: rename tophat output to <sample>.bam (as in the 2018 dirs)
for S in $SAMPLES; do
  [ -f "$S/accepted_hits.bam" ] && mv "$S/accepted_hits.bam" "$S/$S.bam"
done
echo "alignment done."
