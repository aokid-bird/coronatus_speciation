#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def run(cmd, **kwargs):
    print("RUN:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, **kwargs)


def main():
    ap = argparse.ArgumentParser(description="Prepare TreeMix input from ACF using glactools: rename, optional meld, optional root, and acf2treemix")
    ap.add_argument("--acf", required=True, help="Path to input .acf.gz from glactools vcfm2acf")
    ap.add_argument("--bamlist", required=True, help="Bamlist used in ANGSD (lines correspond to samples)")
    ap.add_argument("--samples", required=True, help="TSV with columns including 'sample' and the group column")
    ap.add_argument("--group-col", required=True, help="Column name in samples TSV for population/group")
    ap.add_argument("--mode", choices=["pop", "indiv"], default="pop")
    ap.add_argument("--include-outgroups", action="store_true", help="Include outgroups and set root/anc")
    ap.add_argument("--outgroups-tsv", default=None, help="TSV with outgroup 'sample_id' and 'taxon' columns")
    ap.add_argument("--root-label", default="", help="Root population label (if empty, derive from outgroups.tsv taxon when include-outgroups is set)")
    ap.add_argument("--merge-outgroups", choices=["true", "false"], default="false", help="If 'true', merge multiple outgroup taxa into a single root group; if 'false', keep separate taxa")
    ap.add_argument("--glactools-bin", default="workflow/bin/glactools", help="Path to glactools binary")
    ap.add_argument("--outdir", required=True, help="Output directory for treemix prep")
    ap.add_argument("--treemix-input", required=True, help="Output treemix matrix path (uncompressed)")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Build rename mapping old->new
    bamlist = Path(args.bamlist).read_text().strip().splitlines()
    old_names = [b.strip() for b in bamlist if b.strip()]
    # sample name = bam basename without .bam
    samples_from_bam = [Path(p).stem for p in old_names]

    # Verify ACF sample order via glactools view -p (best effort)
    try:
        res = subprocess.run([args.glactools_bin, "view", "-p", args.acf], check=True, capture_output=True, text=True)
        acf_names = [line.strip() for line in res.stdout.splitlines() if line.strip()]
        # If counts match and the last n lines equal old_names, use those as old_names
        # glactools prints the list; often matches bamlist. If not equal length, fall back.
        if len(acf_names) == len(old_names):
            old_names = acf_names
    except Exception:
        pass

    # Build new names initially equal to sample IDs
    new_names = samples_from_bam.copy()

    # 2) Prepare merge (meld) groups if mode==pop
    samples_df = pd.read_csv(args.samples, sep="\t")
    if "sample" not in samples_df.columns:
        raise SystemExit("samples TSV must include a 'sample' column")
    if args.group_col not in samples_df.columns:
        raise SystemExit(f"samples TSV must include the group column: {args.group_col}")

    group_map = dict(zip(samples_df["sample"].astype(str), samples_df[args.group_col].astype(str)))

    # If including outgroups, add them and map to root label for group
    og_merge_label = None
    og_taxon_labels = {}
    if args.include_outgroups:
        if not args.outgroups_tsv:
            raise SystemExit("--outgroups-tsv is required when --include-outgroups is set")
        og = pd.read_csv(args.outgroups_tsv, sep="\t")
        if "sample_id" not in og.columns:
            raise SystemExit("outgroups TSV must include a 'sample_id' column")

        def abbr_taxon(t):
            if pd.isna(t):
                return None
            words = str(t).strip().split()
            return "".join([w[:2] for w in words if len(w) > 0])

        # Determine labeling for outgroups
        merge_ogs = (args.merge_outgroups.lower() == "true")
        have_taxon = ("taxon" in og.columns)
        if have_taxon:
            og["taxon_label"] = og["taxon"].apply(abbr_taxon)
            # Map sample->taxon_label for later grouping
            og_taxon_labels = dict(zip(og["sample_id"].astype(str), og["taxon_label"].astype(str)))

        if merge_ogs:
            # collapse all outgroups to one label
            if args.root_label:
                og_merge_label = args.root_label
            else:
                # single taxon → use that as label; else fall back to OUT
                if have_taxon:
                    uniq = sorted([x for x in og["taxon_label"].dropna().unique().tolist() if x])
                    og_merge_label = uniq[0] if len(uniq) == 1 else "OUT"
                else:
                    og_merge_label = "OUT"
            for sid in og["sample_id"].astype(str).tolist():
                group_map[str(sid)] = og_merge_label
        else:
            # retain separate taxa (or per-sample if taxon missing)
            for sid in og["sample_id"].astype(str).tolist():
                if have_taxon and og_taxon_labels.get(str(sid)) not in (None, "", "nan"):
                    group_map[str(sid)] = og_taxon_labels[str(sid)]
                else:
                    group_map[str(sid)] = str(sid)

    # 3) Execute glactools pipeline
    tmp1 = outdir / "gl_treemix_tmp1.acf.gz"
    tmp2 = outdir / "gl_treemix_tmp2.acf.gz"
    tmp3 = outdir / "gl_treemix_tmp3.acf.gz"
    final_acf = outdir / "input.acf.gz"

    # rename
    run([args.glactools_bin, "rename", args.acf,
         ",".join(old_names), ",".join(new_names)], stdout=open(tmp1, "wb"))

    # meld if needed
    if args.mode == "pop":
        # Build meld spec: for each group label, join samples that belong to it
        from_to_pairs = []
        # Use only samples that we actually renamed into (new_names)
        new_set = set(new_names)
        groups_to_members = {}
        for s in new_set:
            g = group_map.get(s)
            if g is None:
                continue
            groups_to_members.setdefault(g, []).append(s)
        # build meld command string: "\"a,b\" \"A\" \"c\" \"B\" ..."
        meld_args = []
        for g, members in groups_to_members.items():
            members_csv = ",".join(members)
            meld_args.extend([members_csv, g])
        run([args.glactools_bin, "meld", str(tmp1)] + meld_args, stdout=open(tmp2, "wb"))
    else:
        # indiv mode: if include_outgroups, meld outgroup individuals by taxon when keeping separate taxa,
        # or merge all OG into one label when merge_ogs=true; ingroup remain individuals.
        if args.include_outgroups:
            og_df = pd.read_csv(args.outgroups_tsv, sep="\t")
            og_samples = set(og_df["sample_id"].astype(str).tolist())
            # Build groupings only for outgroup members present in dataset
            new_set = set(new_names)
            present_ogs = [s for s in new_set if s in og_samples]
            if len(present_ogs) == 0:
                tmp2.write_bytes(tmp1.read_bytes())
            else:
                meld_args = []
                if merge_ogs:
                    meld_args.extend([",".join(present_ogs), og_merge_label])
                else:
                    # group by taxon_label (or sample if missing)
                    # rebuild groups_to_members only for OG
                    groups_to_members = {}
                    for s in present_ogs:
                        g = group_map.get(s, s)
                        groups_to_members.setdefault(g, []).append(s)
                    for g, members in groups_to_members.items():
                        members_csv = ",".join(members)
                        meld_args.extend([members_csv, g])
                run([args.glactools_bin, "meld", str(tmp1)] + meld_args, stdout=open(tmp2, "wb"))
        else:
            tmp2.write_bytes(tmp1.read_bytes())

    # set root/anc if including outgroups
    if args.include_outgroups:
        # Determine root label: prefer explicit flag, else derived.
        if args.root_label:
            root_lab = args.root_label
        elif og_merge_label:
            root_lab = og_merge_label
        else:
            # choose the outgroup taxon label with most members (present in dataset)
            og_df = pd.read_csv(args.outgroups_tsv, sep="\t")
            og_samples = og_df["sample_id"].astype(str).tolist()
            # count present members per group label among OG
            counts = {}
            for s in og_samples:
                g = group_map.get(s, None)
                if g is None:
                    continue
                counts[g] = counts.get(g, 0) + 1
            root_lab = max(counts, key=counts.get) if counts else None
        if not root_lab:
            raise SystemExit("Root label could not be determined for outgroups. Provide treemix.root_label or ensure outgroups.tsv has 'taxon'.")
        run([args.glactools_bin, "usepopsrootanc", str(tmp2), root_lab, root_lab], stdout=open(tmp3, "wb"))
        run([args.glactools_bin, "replaceanc", str(tmp2), str(tmp3)], stdout=open(final_acf, "wb"))
    else:
        final_acf.write_bytes(tmp2.read_bytes())

    # acf2treemix (no root, no private)
    with open(args.treemix_input, "wb") as outf:
        run([args.glactools_bin, "acf2treemix", "--noroot", "--noprivate", str(final_acf)], stdout=outf)

    # version check (print to stderr)
    try:
        run([args.glactools_bin, "view", "-P", str(final_acf)])
    except Exception:
        pass


if __name__ == "__main__":
    main()
