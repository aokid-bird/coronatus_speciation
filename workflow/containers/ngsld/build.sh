#!/bin/bash
# Build script for ngsRelate container
set -euo pipefail

# make dir if needed
TARGET_DIR="$HOME/envs/singularity"
IMAGE_NAME="ngsld_1.2.0.sif"
DOCKER_SRC="docker://ghcr.io/zjnolen/ngsld:1.2.0"

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
