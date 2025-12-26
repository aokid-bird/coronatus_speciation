#!/usr/bin/env python3
"""Write SNAPP constraint file from config-defined distributions."""
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--type", required=True, help="Constraint placement (e.g., crown, stem)")
    p.add_argument("--taxa", required=True, help="Comma-separated list of taxa")
    p.add_argument("--distribution", action="append", required=True,
                   help="Distribution string, e.g., lognormal(0,5.11,0.3). Repeatable.")
    p.add_argument("--output", required=True, help="Output constraint filename")
    return p.parse_args()


def main():
    args = parse_args()
    lines = [f"{dist}\t{args.type}\t{args.taxa}" for dist in args.distribution]
    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as handle:
        for line in lines:
            handle.write(f"{line}\n")


if __name__ == "__main__":
    main()
