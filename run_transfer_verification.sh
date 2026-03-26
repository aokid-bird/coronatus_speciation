#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash run_transfer_verification.sh [--report <tsv>] <manifest_tsv> [destination_root]

Behavior:
  - Verifies transfer state using rsync dry-run comparisons.
  - Works with:
      1. data/ssd_mirroring.tsv
         columns: cold_path, ssd_path
      2. transfer/<output_prefix>/cold_storage_paths.tsv
         columns: kind, category, local_path, relative_subpath
         requires destination_root
  - Writes a TSV report summarizing whether each row is in sync.

Examples:
  bash run_transfer_verification.sh data/ssd_mirroring.tsv
  bash run_transfer_verification.sh transfer/defaults/cold_storage_paths.tsv /Volumes/cold_storage/project_a
EOF
}

REPORT_TSV=""
while [[ $# -gt 0 ]]; do
    case "${1}" in
        --report)
            REPORT_TSV="${2:-}"
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

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 1
fi

MANIFEST_TSV="${1}"
DEST_ROOT="${2:-}"

if [[ ! -f "${MANIFEST_TSV}" ]]; then
    echo "Error: manifest does not exist: ${MANIFEST_TSV}" >&2
    exit 1
fi

if [[ -z "${REPORT_TSV}" ]]; then
    REPORT_TSV="$(dirname "${MANIFEST_TSV}")/verification_$(date '+%Y%m%d_%H%M%S').tsv"
fi

header="$(head -n 1 "${MANIFEST_TSV}")"
mode=""
if [[ "${header}" == $'cold_path\tssd_path' ]]; then
    mode="ssd"
elif [[ "${header}" == $'kind\tcategory\tlocal_path\trelative_subpath' ]]; then
    mode="cold"
    if [[ -z "${DEST_ROOT}" ]]; then
        echo "Error: destination_root is required for cold_storage_paths.tsv verification." >&2
        exit 1
    fi
else
    echo "Error: unsupported manifest format: ${MANIFEST_TSV}" >&2
    exit 1
fi

echo -e "source_path\tdestination_path\tstatus\tdetails" > "${REPORT_TSV}"

in_sync_count=0
needs_sync_count=0
missing_count=0
checked_count=0

check_pair() {
    local source_path="$1"
    local destination_path="$2"
    local rsync_output=""
    local status=""
    local details=""

    if [[ ! -e "${source_path}" ]]; then
        status="source_missing"
        details="source_not_found"
        missing_count=$((missing_count + 1))
    elif [[ ! -e "${destination_path}" ]]; then
        status="destination_missing"
        details="destination_not_found"
        missing_count=$((missing_count + 1))
    else
        if [[ -d "${source_path}" ]]; then
            rsync_output="$(rsync -ain --delete --itemize-changes "${source_path}/" "${destination_path}/" 2>&1 || true)"
        else
            rsync_output="$(rsync -ain --itemize-changes "${source_path}" "${destination_path}" 2>&1 || true)"
        fi
        if [[ -z "${rsync_output}" ]]; then
            status="in_sync"
            details="no_changes"
            in_sync_count=$((in_sync_count + 1))
        else
            status="needs_sync"
            details="$(printf '%s' "${rsync_output}" | tr '\n' ';' | sed 's/[[:space:]]\+/ /g' | sed 's/;$/ /')"
            needs_sync_count=$((needs_sync_count + 1))
        fi
    fi

    checked_count=$((checked_count + 1))
    echo -e "${source_path}\t${destination_path}\t${status}\t${details}" >> "${REPORT_TSV}"
}

{
    read -r _
    while IFS=$'\t' read -r c1 c2 c3 c4; do
        [[ -z "${c1}" ]] && continue
        if [[ "${mode}" == "ssd" ]]; then
            check_pair "${c1}" "${c2}"
        else
            check_pair "${c3}" "${DEST_ROOT}/${c4}"
        fi
    done
} < "${MANIFEST_TSV}"

echo "Manifest: ${MANIFEST_TSV}"
echo "Report: ${REPORT_TSV}"
echo "Checked: ${checked_count}"
echo "In sync: ${in_sync_count}"
echo "Needs sync: ${needs_sync_count}"
echo "Missing: ${missing_count}"
