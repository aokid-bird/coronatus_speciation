import argparse
import pandas as pd
import gzip
from functools import reduce

parser = argparse.ArgumentParser(description="***** This script takes multiple ANGSD .geno.gz files, extracts positions,\
                                 and finds intersecting sites. *****")
parser.add_argument("--genofiles", nargs='+', required=True, help="List of geno.gz files")
parser.add_argument("--out_sites", required=True, help="Output file for intersecting sites")
parser.add_argument("--out_scafs", required=True, help="Output file for scaffolds")
args = parser.parse_args()

def read_sites(filepath):
    with gzip.open(filepath, 'rt') as f:
        df = pd.read_csv(f, sep='\t', usecols=[0, 1], header=None, names=["scaf", "pos"])
        df["tag"] = df["scaf"].astype(str) + ":" + df["pos"].astype(str)
        return df[["tag"]]

sites_list = [read_sites(f) for f in args.genofiles]
intersect = reduce(lambda x, y: pd.merge(x, y, on="tag"), sites_list)
intersect[["scaf", "pos"]] = intersect["tag"].str.split(":", expand=True)

intersect.sort_values(["scaf", "pos"], inplace=True)
intersect[["scaf", "pos"]].to_csv(args.out_sites, sep='\t', index=False, header=False)
intersect[["scaf"]].drop_duplicates().to_csv(args.out_scafs, sep='\t', index=False, header=False)
