#!/bin/bash
# Build script for ngsRelate container
set -euo pipefail

# make dir if needed
TARGET_DIR="$HOME/envs/singularity"
IMAGE_NAME="ngsrelate_20220925.sif"
DOCKER_SRC="docker://ghcr.io/zjnolen/ngsrelate:20220925-ec95c8f"

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

# build singularity image
echo "Building Singularity image..."
singularity build "$IMAGE_NAME" "$DOCKER_SRC"

# generate checksum
echo "Generating checksum..."
sha256sum "$IMAGE_NAME" > "${IMAGE_NAME}.sha256"

echo "Done."
echo "Image stored at: $TARGET_DIR/$IMAGE_NAME"
echo "Checksum stored at: $TARGET_DIR/${IMAGE_NAME}.sha256"
