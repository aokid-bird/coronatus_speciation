#!/usr/bin/env python3
"""Extract contig lengths from a FASTA index (.fai) and export as CSV."""

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True, help="Path to the FASTA index (.fai)")
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output CSV file (columns: accession,length)",
    )
    return parser.parse_args()


def read_fai(fai_path: Path):
    with fai_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            accession, length = fields[0], fields[1]
            yield accession, int(length)


def write_csv(rows, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "length"])
        for accession, length in rows:
            writer.writerow([accession, length])


def main() -> None:
    args = parse_args()
    fai_path = Path(args.fai)
    out_path = Path(args.output)
    rows = list(read_fai(fai_path))
    write_csv(rows, out_path)


if __name__ == "__main__":
    main()
