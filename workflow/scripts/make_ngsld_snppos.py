# workflow/scripts/make_ngsld_snppos.py
#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--beagle", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    
    # assuming 1st row for header
    # read as chunks
    reader = pd.read_csv(
        args.beagle,
        sep="\t",
        compression="gzip",
        header=None,
        skiprows=1,
        chunksize=100000,
        usecols=[0]      # read only V1 column
    )

    mode = "w"
    for chunk in reader:
        # "scaf_pos" → ['scaf','snppos']
        df2 = chunk.iloc[:,0] \
            .str.split(pat="_", n=1, expand=True) \
            .rename(columns={0:"scaf", 1:"snppos"})
        df2["snppos"] = df2["snppos"].astype(int)

        # without headers
        # mode = "w" for the first time to write out
        # and then append after mode = "a"
        df2.to_csv(
            args.out,
            sep="\t",
            header=False,
            index=False,
            mode=mode
        )
        mode = "a"

if __name__ == "__main__":
    main()
