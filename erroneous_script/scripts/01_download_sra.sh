#!/bin/bash
# =============================================================================
# 01_download_sra.sh - Download the 6 RNA-seq libraries from ENA/NCBI
# =============================================================================
# Sample manifest (ENA/GEO metadata). NOTE on naming: the GEO submission
# labels these runs "_base_" and that is the tissue used throughout; the
# original 2018 local analysis had the base/seedling file names swapped.
#
#   SRR       GEO sample title                  used here as
#   SRR9070097  WT_base_mRNA180209_rep1       -> WT_base_mRNA180209_rep1
#   SRR9070098  WT_base_mRNA180209_rep2       -> WT_base_mRNA180209_rep2
#   SRR9070099  WT_base_mRNA180209_rep3       -> WT_base_mRNA180209_rep3
#   SRR9070106  osnrpd1ab_base_mRNA180209_rep1 -> pol4ab_base_mRNA180209_rep1
#   SRR9070107  osnrpd1ab_base_mRNA180209_rep2 -> pol4ab_base_mRNA180209_rep2
#   SRR9070108  osnrpd1ab_base_mRNA180209_rep3 -> pol4ab_base_mRNA180209_rep3
#
# The 6 libraries are the ones analysed in the original paper (verified:
# read counts and read sequences of the SRR runs are identical to the
# original local FASTQ files).
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/rawdata" && cd "$ROOT/rawdata"

BASE=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR907
declare -A MAP=(
  [SRR9070097]=WT_base_mRNA180209_rep1
  [SRR9070098]=WT_base_mRNA180209_rep2
  [SRR9070099]=WT_base_mRNA180209_rep3
  [SRR9070106]=pol4ab_base_mRNA180209_rep1
  [SRR9070107]=pol4ab_base_mRNA180209_rep2
  [SRR9070108]=pol4ab_base_mRNA180209_rep3
)
declare -A SUB=( [SRR9070097]=007 [SRR9070098]=008 [SRR9070099]=009
                 [SRR9070106]=006 [SRR9070107]=007 [SRR9070108]=008 )

# ENA fastq_md5 (verified during this reproduction):
declare -A MD5=(
  [SRR9070097]="ea6682133a01fa35448f9d182109ec27;32205cf888beff17fbc0c7e06c8f2e35"
  [SRR9070098]="70283cccaf90d2001d31d6b59c720ce0;08ca56aee0f059cf8bb4d8f7018046eb"
  [SRR9070099]="967678b3525ca0cd5b932861e0aaf227;bda1c72bc946e8eb0caed491a9c913b4"
  [SRR9070106]="c8ed53760dd8b3fe03bd5fdec7de28a7;4a0763340d68bb518809a8c32a1afd7d"
  [SRR9070107]="e3ff23125163a5374d474da40add5a75;db2f0183cdcc6ca539b4bf6d454eb723"
  [SRR9070108]="1422f17817186eed36bf0b1a3747beb5;c7695f417bce041d31053778f3d0752d"
)

for SRR in "${!MAP[@]}"; do
  for MATE in 1 2; do
    URL="$BASE/${SUB[$SRR]}/$SRR/${SRR}_${MATE}.fastq.gz"
    F="${SRR}_${MATE}.fastq.gz"
    [ -f "$F" ] || aria2c -x 8 -s 8 -c --max-tries=10 --retry-wait=15 -o "$F" "$URL"
  done
done

# verify md5 (all 12 files)
for SRR in "${!MAP[@]}"; do
  M=(${MD5[$SRR]//;/ })
  echo "${M[0]}  ${SRR}_1.fastq.gz" >> md5check.txt
  echo "${M[1]}  ${SRR}_2.fastq.gz" >> md5check.txt
done
md5sum -c md5check.txt && rm md5check.txt

# rename to the sample names used by the pipeline (symlinks)
ln -sf SRR9070097_1.fastq.gz WT_base_mRNA180209_rep1_PE1.fq.gz
ln -sf SRR9070097_2.fastq.gz WT_base_mRNA180209_rep1_PE2.fq.gz
ln -sf SRR9070098_1.fastq.gz WT_base_mRNA180209_rep2_PE1.fq.gz
ln -sf SRR9070098_2.fastq.gz WT_base_mRNA180209_rep2_PE2.fq.gz
ln -sf SRR9070099_1.fastq.gz WT_base_mRNA180209_rep3_PE1.fq.gz
ln -sf SRR9070099_2.fastq.gz WT_base_mRNA180209_rep3_PE2.fq.gz
ln -sf SRR9070106_1.fastq.gz pol4ab_base_mRNA180209_rep1_PE1.fq.gz
ln -sf SRR9070106_2.fastq.gz pol4ab_base_mRNA180209_rep1_PE2.fq.gz
ln -sf SRR9070107_1.fastq.gz pol4ab_base_mRNA180209_rep2_PE1.fq.gz
ln -sf SRR9070107_2.fastq.gz pol4ab_base_mRNA180209_rep2_PE2.fq.gz
ln -sf SRR9070108_1.fastq.gz pol4ab_base_mRNA180209_rep3_PE1.fq.gz
ln -sf SRR9070108_2.fastq.gz pol4ab_base_mRNA180209_rep3_PE2.fq.gz
echo "rawdata ready."
