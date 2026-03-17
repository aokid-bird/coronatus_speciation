#!/bin/bash
# run_pipeline.sh

# stop when error comes out
set -euo pipefail

# activate environment
source /app/miniconda/22.11.1-1/python_3.10/etc/profile.d/conda.sh
conda activate bioinfo_pipeline

# absolute paths and variables
CONDA_PREFIX_SNAKEMAKE="$HOME/envs/conda"
SINGULARITY_PREFIX_SNAKEMAKE="$HOME/envs/singularity"
CONFIG=config/config.yaml
PROFILE=profile/default
LOGDIR=$(python workflow/scripts/export_paths.py $CONFIG)
echo "Your LOGDIR is \"$LOGDIR\""
mkdir -p $LOGDIR

# export Slurm profile settings from config
eval "$(
python - "$CONFIG" "$LOGDIR" <<'PY'
import shlex
import sys
import yaml

config_path, logdir = sys.argv[1], sys.argv[2]
with open(config_path) as fh:
    config = yaml.safe_load(fh) or {}

slurm = config.get("slurm", {}) or {}
slurm_logdir = f"{logdir}/slurm"

extra_args = []
for cli_opt, key in (
    ("--account", "account"),
    ("--qos", "qos"),
    ("--reservation", "reservation"),
    ("--mail-type", "mail_type"),
    ("--mail-user", "mail_user"),
):
    value = slurm.get(key)
    if value:
        extra_args.append(f"{cli_opt}={value}")

env = {
    "SNAKEMAKE_SLURM_PARTITION": str(slurm.get("partition", "default")),
    "SNAKEMAKE_SLURM_LOGDIR": slurm_logdir,
    "SNAKEMAKE_SLURM_EXTRA_ARGS": " ".join(extra_args),
}

for key, value in env.items():
    print(f"export {key}={shlex.quote(value)}")
PY
)"
mkdir -p "$SNAKEMAKE_SLURM_LOGDIR"
echo "Using Slurm partition \"$SNAKEMAKE_SLURM_PARTITION\""

# create conda environment
echo "Creating conda environments..."
snakemake \
  --use-conda \
  --conda-create-envs-only \
  --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
  --rerun-incomplete \
  --conda-frontend conda \
  --unlock \
  --cores 1

# build containers
echo "Building containers..."
for script in workflow/containers/*/build.sh; do
    # obtain image name from build.sh
    image_name=$(grep '^IMAGE_NAME=' "$script" | cut -d= -f2- | tr -d '"')
    # if image name does not exist
    if [ -z "$image_name" ]; then
      echo "Warning: IMAGE_NAME not found in $script. Running script anyway."
      bash "$script"
      continue
    fi
    # generate path for the sif image
    sif_path="$SINGULARITY_PREFIX_SNAKEMAKE/$image_name"
    if [ -f "$sif_path" ]; then
        echo "Skipping $script (already built: $sif_path)"
        continue
    fi
    echo "Running $script"
    bash "$script"
done

TARGETS=("$@")

# create DAG (for specified targets if provided)
if [ ${#TARGETS[@]} -gt 0 ]; then
  snakemake --dag --rerun-incomplete --unlock --snakefile workflow/Snakefile --configfile $CONFIG "${TARGETS[@]}" | dot -Tpdf > "${LOGDIR}/dag.pdf"
else
  snakemake --dag --rerun-incomplete --unlock | dot -Tpdf > "${LOGDIR}/dag.pdf"
fi

# Snakemake dryrun
snakemake \
    --dryrun \
    --printshellcmds \
    --cores 6 \
    --snakefile workflow/Snakefile \
    --rerun-incomplete \
    --configfile config/config.yaml \
    --unlock \
    "${TARGETS[@]}" > ${LOGDIR}/dryrun.log 2>&1

# Snakemake run
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --profile "$PROFILE" \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
    --use-singularity \
    --rerun-incomplete \
    "${TARGETS[@]}" # optional targets for partial runs

# For partial debug
# snakemake \
#     --snakefile workflow/Snakefile \
#     --configfile config/config.yaml \
#     -R intersect_sites_group \
#     --profile "$PROFILE" \
#     --use-conda \
#     --conda-frontend conda \
#     --use-singularity
