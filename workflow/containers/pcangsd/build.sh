#!/bin/bash
#workflow/containers/pcangsd/build.sh
set -euo pipefail

TARGET_DIR="$HOME/envs/singularity"
IMAGE_NAME="pcangsd_1.35.sif"
DOCKER_SRC="docker://ghcr.io/zjnolen/pcangsd:1.35"

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

echo "Building Singularity image..."
singularity build "$IMAGE_NAME" "$DOCKER_SRC"

echo "Generating checksum..."
sha256sum "$IMAGE_NAME" > "${IMAGE_NAME}.sha256"

echo "Done."
