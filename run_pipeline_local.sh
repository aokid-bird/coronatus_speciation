#!/bin/bash
# run_pipeline.sh

set -euo pipefail

# activate environment
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    echo "Unable to initialize conda shell integration." >&2
    exit 1
fi
conda activate bioinfo_pipeline

# absolute paths and variables
CONDA_PREFIX_SNAKEMAKE="$HOME/envs/conda"
CONFIG=config/config.yaml
CONDARC_SNAKEMAKE="$(pwd)/conda/condarc_snakemake.yaml"
LOGDIR=$(python workflow/scripts/export_paths.py "$CONFIG")
mkdir -p "$LOGDIR"
echo "Your LOGDIR is \"$LOGDIR\""

# Keep R package resolution inside each conda environment.
unset R_LIBS_USER R_PROFILE_USER R_ENVIRON_USER

# Force Snakemake-created conda envs to ignore user/system default channels.
export CONDARC="$CONDARC_SNAKEMAKE"

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
    --configfile "$CONFIG" > "${LOGDIR}/dryrun.log" 2>&1

# for SRA downlaod/reference download, do the following in the local MacStudio
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --config environment=local \
    --use-conda \
    --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
    --conda-frontend conda \
    --cores 14
