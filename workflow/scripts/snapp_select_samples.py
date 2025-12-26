#!/usr/bin/env python3
"""Select SNAPP individuals with lowest missingness per population."""
import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno", required=True, help="ANGSD genotype file (.geno.gz)")
    p.add_argument("--bamlist", required=True, help="Bamlist used for the ANGSD run")
    p.add_argument("--samples", required=True, help="Sample metadata TSV")
    p.add_argument("--group-col", required=True, help="Column in metadata with population labels")
    p.add_argument("--populations", required=True, nargs="+", help="Target populations to evaluate")
    p.add_argument("--max-per-pop", type=int, default=4, help="Maximum number of individuals per population to keep")
    p.add_argument("--missingness-threshold", type=float, default=None,
                   help="Optional maximum fraction of missing genotypes (0-1). Individuals above are dropped")
    p.add_argument("--exclude-sample", action="append", default=[],
                   help="Sample ID to exclude from consideration (may be provided multiple times)")
    p.add_argument("--summary-tsv", required=True, help="Output TSV summarising missingness per sample")
    p.add_argument("--selected-list", required=True, help="Output text file listing retained sample IDs")
    return p.parse_args()


def _read_bamlist(bamlist_path: Path) -> List[str]:
    names: List[str] = []
    with open(bamlist_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            stem = Path(line).name
            for suffix in (".bam", ".cram"):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
            if stem.endswith(".sam"):
                stem = stem[: -4]
            # remove trailing patterns like _sliceXYZ (optional)
            if "_slice" in stem:
                stem = stem.split("_slice", 1)[0]
            names.append(stem)
    return names


def _load_metadata(samples_path: Path, group_col: str) -> pd.DataFrame:
    df = pd.read_csv(samples_path, sep="\t")
    if group_col not in df.columns:
        raise ValueError(f"Column '{group_col}' not found in {samples_path}")
    return df


def _count_missing(geno_path: Path, analysis_samples: List[str], tracked: List[str]) -> Tuple[Dict[str, int], int]:
    """Return missingness counts for tracked samples and total site count."""
    analysis_index = {sample: idx for idx, sample in enumerate(analysis_samples)}
    for sample in tracked:
        if sample not in analysis_index:
            raise ValueError(
                f"Sample '{sample}' not present in ANGSD bamlist; available: {analysis_samples}"
            )

    missing_counts = {sample: 0 for sample in tracked}
    total_sites = 0
    offset = None
    positions: Dict[str, int] = {}

    with gzip.open(geno_path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            # Drop trailing empty column caused by ending tab if present
            if row and row[-1] == "":
                row = row[:-1]
            if offset is None:
                offset = len(row) - len(analysis_samples)
                if offset < 0:
                    raise ValueError(
                        f"Row has fewer columns ({len(row)}) than bamlist entries ({len(analysis_samples)})."
                    )
                positions = {sample: offset + analysis_index[sample] for sample in tracked}
            total_sites += 1
            for sample, col_idx in positions.items():
                if col_idx >= len(row):
                    raise ValueError(
                        f"Column index {col_idx} for sample {sample} exceeds row length {len(row)}"
                    )
                val = row[col_idx].strip()
                if not val:
                    missing_counts[sample] += 1
                    continue
                upper = val.upper()
                if "N" in upper or upper in {"-1", "NA"}:
                    missing_counts[sample] += 1
    return missing_counts, total_sites


def main() -> None:
    args = parse_args()

    geno_path = Path(args.geno)
    analysis_bamlist = Path(args.bamlist)
    samples_path = Path(args.samples)

    if args.max_per_pop < 1:
        raise ValueError("--max-per-pop must be >= 1")

    analysis_samples = _read_bamlist(analysis_bamlist)
    meta_all = _load_metadata(samples_path, args.group_col)
    exclude = {str(s).strip() for s in args.exclude_sample if str(s).strip()}
    missing_excludes = sorted(exclude - set(meta_all["sample"].astype(str))) if exclude else []
    if missing_excludes:
        print(
            "[snapp_select_samples] Warning: requested exclusions not found in metadata: "
            + ", ".join(missing_excludes),
            file=sys.stderr,
        )

    meta = meta_all[meta_all[args.group_col].isin(args.populations)].copy()
    if exclude:
        meta = meta[~meta["sample"].astype(str).isin(exclude)].copy()
    if meta.empty:
        raise ValueError("No samples remain after filtering metadata by target populations")

    tracked_samples = meta["sample"].astype(str).tolist()
    missing_counts, total_sites = _count_missing(
        geno_path, analysis_samples, tracked_samples
    )
    if total_sites == 0:
        raise ValueError(f"No sites found in genotype file {geno_path}")

    records = []
    for pop in args.populations:
        subset = meta[meta[args.group_col] == pop].copy()
        if subset.empty:
            continue
        subset["missing_count"] = subset["sample"].map(missing_counts)
        subset["missing_prop"] = subset["missing_count"] / float(total_sites)
        subset.sort_values(["missing_prop", "sample"], inplace=True)
        subset["rank"] = range(1, len(subset) + 1)
        subset["selected"] = subset["rank"] <= args.max_per_pop
        if args.missingness_threshold is not None:
            threshold = float(args.missingness_threshold)
            subset.loc[subset["missing_prop"] > threshold, "selected"] = False
        records.append(subset)
    if records:
        result = pd.concat(records, ignore_index=True)
    else:
        result = pd.DataFrame(columns=["sample", args.group_col, "missing_count", "missing_prop", "rank", "selected"])

    # Ensure at least one individual retained per population, optionally relaxing threshold
    shortages = [
        pop
        for pop in args.populations
        if not ((result[args.group_col] == pop) & (result["selected"])).any()
    ]

    if shortages and args.missingness_threshold is not None:
        threshold = float(args.missingness_threshold)
        for pop in shortages:
            pop_mask = result[args.group_col] == pop
            if not pop_mask.any():
                continue
            candidates = result.loc[pop_mask].sort_values(
                ["missing_prop", "rank", "sample"], ascending=[True, True, True]
            )
            top_candidates = candidates.head(args.max_per_pop)
            result.loc[top_candidates.index, "selected"] = True
            missing_values = ", ".join(f"{val:.4f}" for val in top_candidates["missing_prop"])
            print(
                f"[snapp_select_samples] Warning: all individuals for population '{pop}' exceed "
                f"missingness threshold {threshold:.3f}; retaining best {len(top_candidates)} "
                f"sample(s) with missingness {missing_values}.",
                file=sys.stderr,
            )
        # Recompute shortages after relaxing the threshold
        shortages = [
            pop
            for pop in args.populations
            if not ((result[args.group_col] == pop) & (result["selected"])).any()
        ]

    if shortages:
        raise RuntimeError(
            "No individuals retained for populations: " + ", ".join(shortages)
        )

    summary_cols = ["sample", args.group_col, "missing_count", "missing_prop", "rank", "selected"]
    result.sort_values([args.group_col, "rank", "sample"], inplace=True)
    result.to_csv(args.summary_tsv, sep="\t", index=False, columns=summary_cols)

    selected_samples = result.loc[result["selected"], "sample"].tolist()
    Path(args.selected_list).parent.mkdir(parents=True, exist_ok=True)
    with open(args.selected_list, "w") as handle:
        for sample in selected_samples:
            handle.write(f"{sample}\n")


if __name__ == "__main__":
    main()
