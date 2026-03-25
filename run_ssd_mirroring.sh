#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash run_ssd_mirroring.sh [--apply] <cold_root> <ssd_root>

Behavior:
  - By default this runs in rsync dry-run mode.
  - Use --apply to perform the actual copy.
  - The script mirrors the contents of <cold_root>/ into <ssd_root>/.

Examples:
  bash run_ssd_mirroring.sh /Volumes/cold_storage/project_a /Volumes/ssd_storage/project_a
  bash run_ssd_mirroring.sh --apply /Volumes/cold_storage/project_a /Volumes/ssd_storage/project_a
EOF
}

APPLY=false
if [[ "${1:-}" == "--apply" ]]; then
    APPLY=true
    shift
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
fi

COLD_ROOT="${1}"
SSD_ROOT="${2}"

if [[ ! -d "${COLD_ROOT}" ]]; then
    echo "Error: cold storage directory does not exist: ${COLD_ROOT}" >&2
    exit 1
fi

mkdir -p "${SSD_ROOT}"

RSYNC_ARGS=(
    -avh
    --info=progress2
    --human-readable
)

if [[ "${APPLY}" == false ]]; then
    RSYNC_ARGS+=(--dry-run)
    echo "Preview mode only. Re-run with --apply to perform the mirroring."
fi

echo "Cold storage: ${COLD_ROOT}"
echo "SSD mirror:   ${SSD_ROOT}"
echo "Command: rsync ${RSYNC_ARGS[*]} ${COLD_ROOT}/ ${SSD_ROOT}/"

rsync "${RSYNC_ARGS[@]}" "${COLD_ROOT}/" "${SSD_ROOT}/"
