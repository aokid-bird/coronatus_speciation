
# scripts/process_unlinked_and_summary.py
#!/usr/bin/env python3
import argparse, subprocess
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--beagle",      required=True)
    p.add_argument("--unlinkedid",  required=True)
    p.add_argument("--snppos",      required=True)
    p.add_argument("--unlinkedbeagle", required=True)
    p.add_argument("--summary",     required=True)
    args = p.parse_args()

    # 1) total number of sites
    with open(args.snppos) as f:
        num_sites = sum(1 for _ in f)
    # 2) load unlinked IDs
    ul = pd.read_csv(args.unlinkedid, header=None)[0].str.replace(":", "_")
    num_sites_unlink = len(ul)

    # 3) read Beagle file as chunks and extract only unlinked site
    reader = pd.read_csv(
        args.beagle, sep="\t", compression="gzip",
        header=None, chunksize=100000
    )
    header = next(reader)
    frames = [header]
    for chunk in reader:
        filtered = chunk[chunk.iloc[:,0].isin(ul)]
        frames.append(filtered)
    result = pd.concat(frames, ignore_index=True)

    # write as a plane text and compress it by gzip
    tmp = args.unlinkedbeagle.rstrip(".gz")
    result.to_csv(tmp, sep="\t", header=False, index=False)
    subprocess.run(["gzip", "-f", tmp], check=True)

    # 4) output summary csv
    pd.DataFrame({
        "num_sites":        [num_sites],
        "num_sites_unlink": [num_sites_unlink]
    }).to_csv(args.summary, index=False)

if __name__ == "__main__":
    main()
