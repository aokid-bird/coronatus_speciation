import argparse
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description="Build union of sites from multiple ANGSD .geno.gz files")
parser.add_argument("--genofiles", nargs='+', required=True, help="List of geno.gz files")
parser.add_argument("--out_sites", required=True, help="Output TSV for union sites (scaf\tpos)")
parser.add_argument("--out_scafs", required=True, help="Output TXT of unique scaffolds (one per line)")
args = parser.parse_args()

def read_tags(fp):
    with gzip.open(fp, 'rt') as f:
        df = pd.read_csv(f, sep='\t', usecols=[0, 1], header=None, names=["scaf", "pos"], dtype={0:str,1:int})
    df["tag"] = df["scaf"].astype(str) + ":" + df["pos"].astype(str)
    return df[["tag"]]

dfs = [read_tags(f) for f in args.genofiles]
union = pd.concat(dfs, ignore_index=True).drop_duplicates()
union[["scaf", "pos"]] = union["tag"].str.split(":", expand=True)
union = union.drop(columns=["tag"])\
             .assign(pos=lambda d: d["pos"].astype(int))\
             .sort_values(["scaf", "pos"])\
             .reset_index(drop=True)

union[["scaf", "pos"]].to_csv(args.out_sites, sep='\t', index=False, header=False)
union[["scaf"]].drop_duplicates().to_csv(args.out_scafs, sep='\t', index=False, header=False)

