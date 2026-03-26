#!/usr/bin/env python3

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


def ensure_parent(path_str: str) -> None:
    Path(path_str).parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("selection_file")
    parser.add_argument("sites_file")
    parser.add_argument("out_csv")
    parser.add_argument("out_pdf")
    args = parser.parse_args()

    ensure_parent(args.out_csv)
    ensure_parent(args.out_pdf)

    selection = np.load(args.selection_file)
    selection = np.asarray(selection).reshape(-1)
    pvals = chi2.sf(selection, df=1)

    sites = [line.strip() for line in Path(args.sites_file).read_text().splitlines() if line.strip()]

    n = min(len(sites), len(pvals))
    selection_df = pd.DataFrame(
        {
            "site": pd.to_numeric(sites[:n], errors="coerce"),
            "p": pvals[:n],
        }
    )
    selection_df["selection"] = (selection_df["p"] < (0.05 / max(n, 1))).astype(int)
    selection_df.to_csv(args.out_csv, index=False)

    if n == 0:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.text(0.5, 0.5, "No selection statistics available", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(args.out_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    observed = np.sort(selection[:n])
    expected_probs = (np.arange(1, n + 1) - 0.5) / n
    expected = chi2.ppf(expected_probs, df=1)
    lambda_gc = np.median(observed) / chi2.ppf(0.5, df=1)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=10)
    upper = max(np.nanmax(expected), np.nanmax(observed))
    ax.plot([0, upper], [0, upper], color="red", linewidth=1.5)
    ax.set_xlabel("Expected")
    ax.set_ylabel("Observed")
    ax.legend([f"lambda={lambda_gc:.2f}"], loc="upper left", frameon=False)
    fig.savefig(args.out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
