#!/usr/bin/env python3
"""Prepare SNAPP sample renaming files from selected individuals and bamlist."""
import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bamlist", required=True, help="SNAPP bamlist with ingroup and outgroup BAMs")
    p.add_argument("--summary", required=True, help="Selected sample summary TSV")
    p.add_argument("--outgroups", required=True, help="Outgroup metadata TSV")
    p.add_argument("--group-col", required=True, help="Column with ingroup population labels")
    p.add_argument("--sample-names-out", required=True, help="Output sample names text file for bcftools reheader")
    p.add_argument("--pop-table", required=True, help="Output TSV mapping species to individuals")
    p.add_argument("--metadata-out", required=True, help="Output TSV with detailed mapping")
    p.add_argument("--population-alias", action="append", default=[],
                   help="Override population alias using POP=ALIAS format")
    p.add_argument("--outgroup-alias", action="append", default=[],
                   help="Override outgroup alias using SAMPLE_ID=ALIAS format")
    return p.parse_args()


def parse_alias_list(entries: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"Alias entry '{entry}' must be in KEY=VALUE format")
        key, value = entry.split("=", 1)
        mapping[key.strip()] = value.strip()
    return mapping


def sanitize(label: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z]+", "_", label)
    cleaned = cleaned.strip("_")
    return cleaned or "X"


def load_outgroup_labels(path: Path, overrides: Dict[str, str]) -> Dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    labels: Dict[str, str] = {}
    for _, row in df.iterrows():
        sample_id = str(row.get("sample_id"))
        if sample_id in overrides:
            labels[sample_id] = sanitize(overrides[sample_id])
        else:
            taxon = str(row.get("taxon", sample_id))
            if taxon and taxon.lower() != "nan":
                labels[sample_id] = sanitize(taxon.replace(" ", ""))
            else:
                labels[sample_id] = sanitize(sample_id)
    return labels


def read_bam_samples(path: Path) -> List[str]:
    samples: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            name = Path(line).name
            if name.endswith(".bam"):
                name = name[:-4]
            if name.endswith(".cram"):
                name = name[:-5]
            if "_slice" in name:
                name = name.split("_slice", 1)[0]
            samples.append(name)
    return samples


def main() -> None:
    args = parse_args()
    summary = pd.read_csv(args.summary, sep="\t")
    summary = summary[summary["selected"]]
    ingroup_map = summary.set_index("sample")[args.group_col].to_dict()
    population_alias = {k: sanitize(v) for k, v in parse_alias_list(args.population_alias).items()}
    outgroup_alias = {k: sanitize(v) for k, v in parse_alias_list(args.outgroup_alias).items()}
    outgroup_labels = load_outgroup_labels(Path(args.outgroups), outgroup_alias)

    pop_to_alias: Dict[str, str] = {}
    for pop in set(ingroup_map.values()):
        if pop in population_alias:
            pop_to_alias[pop] = population_alias[pop]
        else:
            pop_to_alias[pop] = sanitize(str(pop))

    bam_samples = read_bam_samples(Path(args.bamlist))
    counts: Dict[str, int] = {}
    records = []

    for sample in bam_samples:
        if sample in ingroup_map:
            pop = ingroup_map[sample]
            alias = pop_to_alias[pop]
            source = "ingroup"
        elif sample in outgroup_labels:
            pop = outgroup_labels[sample]
            alias = pop
            source = "outgroup"
        else:
            raise ValueError(
                f"Sample '{sample}' from bamlist not present in selected ingroup list or outgroup metadata"
            )
        counts.setdefault(alias, 0)
        counts[alias] += 1
        new_name = f"{alias}_{counts[alias]}"
        records.append({
            "original_sample": sample,
            "population": pop,
            "alias": alias,
            "individual": new_name,
            "source": source,
        })

    df = pd.DataFrame.from_records(records)

    sample_names = df["individual"].tolist()
    Path(args.sample_names_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.sample_names_out, "w") as sample_file:
        for name in sample_names:
            sample_file.write(f"{name}\n")

    pop_df = df[["population", "individual"]].copy()
    pop_df.columns = ["species", "individual"]
    pop_df.to_csv(args.pop_table, sep="\t", index=False)

    df.to_csv(args.metadata_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
