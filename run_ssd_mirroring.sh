#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash run_ssd_mirroring.sh [--apply] [--log-dir <dir>] <manifest_tsv>

Behavior:
  - By default this runs in rsync dry-run mode.
  - Use --apply to perform the actual copy.
  - The manifest should usually live at data/ssd_mirroring.tsv.
  - The TSV must contain a header with the columns:
      cold_path    ssd_path
  - Each row mirrors one cold-storage path into the corresponding SSD path.
  - Directory and file behavior are inferred from the cold_path on disk.

Examples:
  bash run_ssd_mirroring.sh data/ssd_mirroring.tsv
  bash run_ssd_mirroring.sh --apply data/ssd_mirroring.tsv
EOF
}

APPLY=false
LOG_DIR=""

while [[ $# -gt 0 ]]; do
    case "${1}" in
        --apply)
            APPLY=true
            shift
            ;;
        --log-dir)
            LOG_DIR="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 1
fi

MANIFEST_TSV="${1}"
if [[ -z "${LOG_DIR}" ]]; then
    LOG_DIR="$(dirname "${MANIFEST_TSV}")/logs"
fi

if [[ ! -f "${MANIFEST_TSV}" ]]; then
    echo "Error: manifest does not exist: ${MANIFEST_TSV}" >&2
    exit 1
fi

mkdir -p "${LOG_DIR}"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="${LOG_DIR}/ssd_mirroring_${TIMESTAMP}.log"

RSYNC_ARGS=(
    -avh
    --human-readable
    --itemize-changes
    --stats
)

if rsync --version 2>/dev/null | head -n 1 | grep -Eq 'version 3\.'; then
    RSYNC_ARGS+=(--info=progress2)
else
    RSYNC_ARGS+=(--progress)
fi

if [[ "${APPLY}" == false ]]; then
    RSYNC_ARGS+=(--dry-run)
    echo "Preview mode only. Re-run with --apply to perform the mirroring."
fi

log() {
    echo "$*" | tee -a "${LOG_FILE}"
}

log "Manifest: ${MANIFEST_TSV}"
log "Log file: ${LOG_FILE}"
log "Mode: $([[ "${APPLY}" == true ]] && echo apply || echo preview)"

{
    read -r _
    while IFS=$'\t' read -r cold_path ssd_path; do
        [[ -z "${cold_path}" ]] && continue
        [[ "${cold_path}" == \#* ]] && continue

        if [[ -d "${cold_path}" ]]; then
            log "Mirroring directory: ${cold_path} -> ${ssd_path}"
            if [[ "${APPLY}" == true ]]; then
                mkdir -p "${ssd_path}"
                rsync "${RSYNC_ARGS[@]}" "${cold_path}/" "${ssd_path}/" >> "${LOG_FILE}" 2>&1
            else
                log "Preview: rsync ${RSYNC_ARGS[*]} ${cold_path}/ ${ssd_path}/"
                rsync "${RSYNC_ARGS[@]}" "${cold_path}/" "${ssd_path}/" >> "${LOG_FILE}" 2>&1
            fi
        elif [[ -f "${cold_path}" ]]; then
            log "Mirroring file: ${cold_path} -> ${ssd_path}"
            if [[ "${APPLY}" == true ]]; then
                mkdir -p "$(dirname "${ssd_path}")"
                rsync "${RSYNC_ARGS[@]}" "${cold_path}" "${ssd_path}" >> "${LOG_FILE}" 2>&1
            else
                log "Preview: rsync ${RSYNC_ARGS[*]} ${cold_path} ${ssd_path}"
                rsync "${RSYNC_ARGS[@]}" "${cold_path}" "${ssd_path}" >> "${LOG_FILE}" 2>&1
            fi
        else
            log "Skipping missing source path: ${cold_path}"
        fi
    done
} < "${MANIFEST_TSV}"
