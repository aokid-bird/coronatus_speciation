#!/bin/bash
# run_pipeline.sh

# stop when error comes out
set -euo pipefail

# activate environment
source /Users/daisukeaoki/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo_pipeline

# absolute paths and variables
CONDA_PREFIX_SNAKEMAKE="$HOME/envs/conda"
CONFIG=config/config.yaml
LOGDIR=$(python workflow/scripts/export_paths.py $CONFIG)
mkdir -p $LOGDIR
echo "Your LOGDIR is \"$LOGDIR\""

# create conda environment
echo "Creating conda environments..."
snakemake \
  --config environment=local\
  --use-conda \
  --conda-create-envs-only \
  --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
  --conda-frontend conda \
  --cores 1

# create DAG
snakemake --config environment=local --dag | dot -Tpdf > "${LOGDIR}/dag.pdf"

snakemake \
    --config environment=local\
    --dryrun \
    --printshellcmds \
    --cores 6 \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml > ${LOGDIR}/dryrun.log 2>&1

# for SRA downlaod/reference download, do the following in the local MacStudio
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --config environment=local \
    --use-conda \
    --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
    --conda-frontend conda \
    --cores 14
