#!/bin/bash
# =============================================================================
# 00_setup_envs.sh - Create the pinned software environments (conda)
# =============================================================================
# Version choices follow the 2018 original analysis as closely as possible:
#   - tophat 2.1.1  (exact version of the original run)
#   - bowtie 1.1.2  (2017-era; the exact 2018 build is not recorded)
#   - subread/featureCounts 1.6.0 (2017-era)
#   - DESeq 1.39.0 (last installable version of the legacy DESeq package)
#
# IMPORTANT compatibility notes (discovered during reproduction):
#   * bowtie >= 1.2.3 rejects empty read inputs (/dev/null), which breaks
#     TopHat's internal SAM-header generation -> use bowtie 1.1.2.
#   * featureCounts >= 2.0.1 rejects the original parameter combination
#     "-P without -B" -> use subread 1.6.0 (or 1.6.2; both give identical
#     counts on the same BAM files).
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$ROOT/envs"

conda create -y -p "$ROOT/envs/tophat"  --override-channels -c conda-forge -c bioconda \
    python=2.7 tophat=2.1.1 bowtie=1.1.2
conda create -y -p "$ROOT/envs/subread" --override-channels -c conda-forge -c bioconda \
    subread=1.6.0
conda create -y -p "$ROOT/envs/deseq"   --override-channels -c conda-forge -c bioconda \
    bioconductor-deseq=1.39.0

"$ROOT/envs/tophat/bin/tophat"  --version
"$ROOT/envs/tophat/bin/bowtie"  --version | head -1
"$ROOT/envs/subread/bin/featureCounts" -v 2>&1 | head -1
"$ROOT/envs/deseq/bin/Rscript" -e 'cat("DESeq", as.character(packageVersion("DESeq")), "\n")'
echo "Environments ready."
