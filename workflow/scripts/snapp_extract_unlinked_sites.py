#!/usr/bin/env python3
"""Extract unlinked SNP coordinates for SNAPP from ngsLD IDs and Beagle markers."""
import argparse
import gzip
from pathlib import Path
from typing import List, Tuple


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beagle", required=True, help="ANGSD Beagle file for SNAPP run")
    p.add_argument("--unlinked-id", required=True, help="List of unlinked site IDs (chrom:pos)")
    p.add_argument("--out-sites", required=True, help="Output TSV with chrom and pos of retained markers")
    p.add_argument("--out-summary", required=True, help="Output TSV summarising counts")
    return p.parse_args()


def load_unlinked_ids(path: Path) -> List[str]:
    ids: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            ids.append(line)
    return ids


def marker_to_coord(marker: str) -> Tuple[str, str]:
    if "_" not in marker:
        raise ValueError(f"Marker '{marker}' does not contain an underscore separator")
    chrom, pos = marker.rsplit("_", 1)
    return chrom, pos


def main() -> None:
    args = parse_args()
    beagle_path = Path(args.beagle)
    unlinked_ids = load_unlinked_ids(Path(args.unlinked_id))
    if not unlinked_ids:
        raise ValueError(f"No lines found in {args.unlinked_id}")

    target_markers = {marker.replace(":", "_") for marker in unlinked_ids}
    coords: List[Tuple[str, str]] = []
    total_markers = 0

    with gzip.open(beagle_path, "rt") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty Beagle file: {beagle_path}")
        for line in handle:
            line = line.strip()
            if not line:
                continue
            total_markers += 1
            marker = line.split("\t", 1)[0]
            if marker in target_markers:
                chrom, pos = marker_to_coord(marker)
                coords.append((chrom, pos))

    target_count = len(coords)
    with open(args.out_sites, "w") as out_handle:
        for chrom, pos in coords:
            out_handle.write(f"{chrom}\t{pos}\n")

    with open(args.out_summary, "w") as summary_handle:
        summary_handle.write("total_markers\tmatching_markers\n")
        summary_handle.write(f"{total_markers}\t{target_count}\n")


if __name__ == "__main__":
    main()
