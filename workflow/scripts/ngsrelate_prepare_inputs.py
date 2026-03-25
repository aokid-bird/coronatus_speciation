#!/usr/bin/env python3

import pandas as pd
import gzip
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bamlist", required=True)
parser.add_argument("--mafs", required=True)
parser.add_argument("--id_out", required=True)
parser.add_argument("--freq_out", required=True)
parser.add_argument("--bam_dir", required=True)
args = parser.parse_args()

# generate id
df = pd.read_csv(args.bamlist, sep=r"\s+", header=None, engine="python")
df.columns = ["path"]
ids = df["path"].str.replace(f"{args.bam_dir}/", "", regex=False).str.replace(".bam", "", regex=False)
Path(args.id_out).write_text("\n".join(ids) + "\n")

# generate Freq file
with gzip.open(args.mafs, "rt") as fin, open(args.freq_out, "w") as fout:
    next(fin)  # skip header
    for line in fin:
        cols = line.strip().split()
        if len(cols) >= 7:
            fout.write(cols[6] + "\n")
