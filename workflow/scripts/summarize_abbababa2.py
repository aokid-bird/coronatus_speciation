import math
import pandas as pd


def _significance_label(pvalue: float) -> str:
    if pd.isna(pvalue):
        return ""
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return ""


def main(snakemake):
    obs_path = snakemake.input["observed"]
    trans_path = snakemake.input["transrem"]
    obs = pd.read_csv(obs_path, sep="\t")
    trans = pd.read_csv(trans_path, sep="\t")

    if obs.empty:
        raise ValueError(f"Observed D-statistics file has no rows: {obs_path}")

    n_pairs = len(obs)
    obs = obs.assign(file="Uncorrected")
    trans = trans.assign(file="Uncorrected Transition Removed")
    df = pd.concat([obs, trans], ignore_index=True)

    var_col = "V(JK-D)"
    if var_col not in df.columns:
        raise ValueError(f"Expected column '{var_col}' not found in D-statistics tables: {obs_path}")

    df["sd"] = df[var_col].clip(lower=0).apply(math.sqrt)
    df["adjsd"] = 3 * df["sd"]
    df["adjpvalue"] = df["pvalue"] * n_pairs
    df["label"] = df["adjpvalue"].apply(_significance_label)
    df["pair"] = df.apply(lambda r: f"D{{{r['H1']},{r['H2']}:{r['H3']}}}", axis=1)

    df.to_csv(snakemake.output["summary"], index=False)


if __name__ == "__main__":
    main(snakemake)
