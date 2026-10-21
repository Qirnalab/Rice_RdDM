#!/bin/bash
# =============================================================================
# 02_prepare_reference.sh - Prepare the reference genome and bowtie1 index
# =============================================================================
# The 2018 original used the index "/Share/home/qiyijun/DBs/rice/ricev7",
# which no longer exists. Its exact reference content was recovered from the
# @SQ header of the original BAM files:
#     13 sequences = chloroplast ChrC (134,525 bp) + chromosomes Chr1..Chr12
#     (no ChrUn / ChrSy). Sequence content is identical to the MSU_osa1r7
#     (RAGP v7 / IRGSP-1.0) assembly.
#
# Provide the genome fasta as ref/ricev7.fa yourself (choose one):
#   a) copy the lab file  /qiyijunlab/02.DB/ricev7/ricev7.fa
#      (md5: 6233113c698b58c7268849d11aabb14b, 357 MB)  -- byte-exact;
#   b) reconstruct from public sources: take ChrC + Chr1..Chr12 from the
#      RAGP v7 assembly (osa1_r7.asm.fa) in that order, one sequence each.
# The bowtie index is then rebuilt with bowtie-build 1.1.2 (era-correct
# version). Index build version does not change alignment results as long
# as the sequence content is identical.
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/ref/ricev7" "$ROOT/alignment"

if [ ! -s "$ROOT/ref/ricev7/ricev7.fa" ]; then
  echo "ERROR: place the reference genome at $ROOT/ref/ricev7/ricev7.fa first (see header of this script)."
  exit 1
fi
grep -c ">" "$ROOT/ref/ricev7/ricev7.fa"   # expect 13

export PATH="$ROOT/envs/tophat/bin:$PATH"
bowtie-build --version | head -1
# NOTE: bowtie-build 1.1.2 has no --threads option (single-threaded, ~20 min)
bowtie-build "$ROOT/ref/ricev7/ricev7.fa" "$ROOT/ref/ricev7/ricev7"

# annotation (custom lab GTF, gene models LOC_Os*, provided in annotation/)
cp "$ROOT/annotation/osativa7_newpol4ab.gtf" "$ROOT/alignment/"
echo "reference ready."
