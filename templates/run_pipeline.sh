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
LOGDIR=$(python workflow/scripts/export_paths.py $CONFIG)
echo "Your LOGDIR is \"$LOGDIR\""
mkdir -p $LOGDIR

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
    --cluster-config config/cluster.json \
    --cluster "/opt/nec/nqsv/bin/qsub -q {cluster.queue} \
      -l cpunum_job={cluster.ppn} \
      -l memsz_job={cluster.mem} \
      -N {cluster.jobname} \
      -o ${LOGDIR}/{rule}.o \
      -e ${LOGDIR}/{rule}.e" \
    --jobs 20 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "$CONDA_PREFIX_SNAKEMAKE" \
    --use-singularity \
    --latency-wait 30 \
    --rerun-incomplete \
    --cluster-cancel "qdel" \
    "${TARGETS[@]}" # optional targets for partial runs

# For partial debug
# snakemake \
#     --snakefile workflow/Snakefile \
#     --configfile config/config.yaml \
#     -R intersect_sites_group \
#     --cluster-config config/cluster.json \
#     --cluster "/opt/nec/nqsv/bin/qsub -q {cluster.queue} \
#       -l cpunum_job={cluster.ppn} \
#       -l memsz_job={cluster.mem} \
#       -N {cluster.jobname} \
#       -o {cluster.logdir}/{rule}.o \
#       -e {cluster.logdir}/{rule}.e" \
#     --jobs 6 \
#     --use-conda \
#     --conda-frontend conda \
#     --use-singularity
