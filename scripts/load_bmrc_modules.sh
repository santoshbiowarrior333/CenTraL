#!/usr/bin/env bash
# Load all CenTraL dependencies via Lmod modules.
#
# This file is preconfigured for the BMRC cluster (Oxford). The module
# names and versions below are what BMRC ships as of 2026. On other
# clusters the names will differ - edit the lines below to match what
# 'module avail' shows on your system.
#
# Source this file from your interactive shell or your SLURM job script
# BEFORE running any CenTraL script:
#
#     source scripts/load_bmrc_modules.sh
#     ./scripts/run_dcs_workflow.sh /path/to/bam_pass
#
# CenTraL scripts themselves do not call module load; they assume the
# tools are already on PATH.

module load SAMtools
module load deepTools
module load dorado
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2

# Quick sanity check so a missing module is obvious right away.
missing=0
for t in samtools bamCoverage dorado Rscript; do
    if ! command -v "$t" >/dev/null 2>&1; then
        echo "WARNING: $t not on PATH after module load" >&2
        missing=$((missing + 1))
    fi
done

if [[ $missing -eq 0 ]]; then
    echo "BMRC modules loaded OK:"
    echo "  samtools    : $(samtools --version | head -1)"
    echo "  bamCoverage : $(bamCoverage --version 2>&1 | head -1)"
    echo "  dorado      : $(dorado --version 2>&1 | head -1)"
    echo "  R           : $(R --version 2>&1 | head -1)"
else
    echo "$missing tool(s) missing - either the module name changed or you are not on BMRC." >&2
    echo "Run 'module avail' and update scripts/load_bmrc_modules.sh accordingly." >&2
fi
