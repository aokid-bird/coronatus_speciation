"""
Shared sample and outgroup metadata helpers.

This file owns:
- ingroup FASTQ metadata resolution
- reusable BAM-list selection helpers
- outgroup metadata parsing and outgroup-selection helpers
"""

_SAMPLES_DF = pd.read_csv(config["samples"], sep="\t")
INGROUP_SAMPLE_IDS = _SAMPLES_DF["sample"].astype(str).tolist()

READS_CFG = config.get("reads", {}) or {}
INGROUP_READS_DIR = resolve_ingroup_reads_dir()
INGROUP_READS_META_CFG = READS_CFG.get("ingroup_metadata", {}) or {}

INGROUP_FASTQ_BASENAME_COL = str(INGROUP_READS_META_CFG.get("basename_col", "sample"))
INGROUP_FASTQ_DIR_COL = str(INGROUP_READS_META_CFG.get("dir_col", "fastq_dir"))
INGROUP_FASTQ_PREFIX_COL = str(INGROUP_READS_META_CFG.get("prefix_col", "fastq_prefix"))
INGROUP_FASTQ_R1_SUFFIX_COL = str(INGROUP_READS_META_CFG.get("r1_suffix_col", "fastq_r1_suffix"))
INGROUP_FASTQ_R2_SUFFIX_COL = str(INGROUP_READS_META_CFG.get("r2_suffix_col", "fastq_r2_suffix"))
INGROUP_FASTQ_EXT_COL = str(INGROUP_READS_META_CFG.get("extension_col", "fastq_extension"))
INGROUP_FASTQ_DEFAULT_PREFIX = str(INGROUP_READS_META_CFG.get("default_prefix", ""))
INGROUP_FASTQ_DEFAULT_R1_SUFFIX = str(INGROUP_READS_META_CFG.get("default_r1_suffix", "_1"))
INGROUP_FASTQ_DEFAULT_R2_SUFFIX = str(INGROUP_READS_META_CFG.get("default_r2_suffix", "_2"))
INGROUP_FASTQ_DEFAULT_EXT = str(INGROUP_READS_META_CFG.get("default_extension", ".fastq.gz"))


# Metadata access helpers for sample TSV rows.
def _sample_value(row, colname, default=""):
    if colname in row.index and pd.notna(row[colname]):
        return str(row[colname])
    return default


def _normalize_fastq_ext(ext):
    value = str(ext or "").strip()
    if not value:
        return ""
    return value if value.startswith(".") else f".{value}"


if len(INGROUP_SAMPLE_IDS) != len(set(INGROUP_SAMPLE_IDS)):
    raise ValueError("The 'sample' column in samples.tsv must contain unique analysis sample IDs.")


# Resolve ingroup FASTQ locations once so rules can stay simple.
_INGROUP_FASTQ_RECORDS = {}
INGROUP_STORAGE_SAMPLE_IDS = []
_INGROUP_ANALYSIS_TO_STORAGE = {}
_INGROUP_STORAGE_TO_ANALYSIS = {}
for _, row in _SAMPLES_DF.iterrows():
    sample_id = str(row["sample"]).strip()
    storage_id = _sample_value(row, INGROUP_FASTQ_BASENAME_COL, sample_id).strip()
    if not storage_id:
        raise ValueError(
            f"Sample '{sample_id}' has an empty preprocessing basename after config/metadata resolution."
        )
    if storage_id in _INGROUP_STORAGE_TO_ANALYSIS:
        raise ValueError(
            f"Preprocessing basename '{storage_id}' is reused by both "
            f"'{_INGROUP_STORAGE_TO_ANALYSIS[storage_id]}' and '{sample_id}'."
        )
    sample_dir = _sample_value(row, INGROUP_FASTQ_DIR_COL, INGROUP_READS_DIR).strip()
    if not sample_dir:
        raise ValueError(
            f"Sample '{sample_id}' has an empty FASTQ directory after config/metadata resolution."
        )
    _INGROUP_ANALYSIS_TO_STORAGE[sample_id] = storage_id
    _INGROUP_STORAGE_TO_ANALYSIS[storage_id] = sample_id
    INGROUP_STORAGE_SAMPLE_IDS.append(storage_id)
    _INGROUP_FASTQ_RECORDS[storage_id] = {
        "dir": sample_dir,
        "basename": storage_id,
        "prefix": _sample_value(row, INGROUP_FASTQ_PREFIX_COL, INGROUP_FASTQ_DEFAULT_PREFIX),
        "r1_suffix": _sample_value(row, INGROUP_FASTQ_R1_SUFFIX_COL, INGROUP_FASTQ_DEFAULT_R1_SUFFIX),
        "r2_suffix": _sample_value(row, INGROUP_FASTQ_R2_SUFFIX_COL, INGROUP_FASTQ_DEFAULT_R2_SUFFIX),
        "extension": _normalize_fastq_ext(
            _sample_value(row, INGROUP_FASTQ_EXT_COL, INGROUP_FASTQ_DEFAULT_EXT)
        ),
    }


def ingroup_storage_sample_id(sample_id):
    sid = str(sample_id)
    if sid not in _INGROUP_ANALYSIS_TO_STORAGE:
        raise KeyError(f"Unknown ingroup analysis sample_id '{sid}'")
    return _INGROUP_ANALYSIS_TO_STORAGE[sid]


def ingroup_analysis_sample_id(storage_id):
    sid = str(storage_id)
    if sid not in _INGROUP_STORAGE_TO_ANALYSIS:
        raise KeyError(f"Unknown ingroup preprocessing basename '{sid}'")
    return _INGROUP_STORAGE_TO_ANALYSIS[sid]


# Return the configured FASTQ path for an ingroup sample/read pair.
def ingroup_fastq_path(storage_id, read):
    sid = str(storage_id)
    if sid not in _INGROUP_FASTQ_RECORDS:
        raise KeyError(f"Unknown ingroup preprocessing basename '{sid}'")
    if str(read) not in {"1", "2"}:
        raise ValueError(f"read must be '1' or '2', got '{read}'")
    record = _INGROUP_FASTQ_RECORDS[sid]
    suffix = record["r1_suffix"] if str(read) == "1" else record["r2_suffix"]
    filename = f"{record['prefix']}{record['basename']}{suffix}{record['extension']}"
    return str(Path(record["dir"]) / filename)


def ingroup_storage_bam_path(storage_id, bam_dir=None):
    target_dir = bam_dir or INGROUP_BAM_STORAGE_DIR
    return f"{target_dir}/{storage_id}.bam"


# Derive a sample ID from a BAM/CRAM/SAM path used in bamlists.
def bam_sample_id(path):
    name = Path(str(path)).name
    for suffix in (".bam", ".cram", ".sam"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if "_slice" in name:
        name = name.split("_slice", 1)[0]
    return name


# Reusable ingroup selection helper for analyses that subset samples.
def ingroup_sample_ids(populations=None, exclude_samples=None):
    df = _SAMPLES_DF
    if populations:
        df = df[df[group_col].isin(populations)]
    excluded = set(_parse_list(exclude_samples))
    sample_ids = df["sample"].astype(str).tolist()
    if excluded:
        sample_ids = [sample_id for sample_id in sample_ids if sample_id not in excluded]
    return sample_ids


def ingroup_bam_paths(sample_ids=None, populations=None, exclude_samples=None, bam_dir=None):
    selected = list(sample_ids) if sample_ids is not None else ingroup_sample_ids(
        populations=populations,
        exclude_samples=exclude_samples,
    )
    target_dir = bam_dir or config_bam_dir
    return [f"{target_dir}/{sample_id}.bam" for sample_id in selected]


# Filter an existing bamlist by sample IDs while preserving original order.
def filter_bam_paths_by_sample_ids(bam_paths, exclude_samples=None, include_samples=None):
    excluded = set(_parse_list(exclude_samples))
    included = None if include_samples is None else {str(sample_id) for sample_id in include_samples}
    filtered = []
    for path in bam_paths:
        sample_id = bam_sample_id(path)
        if included is not None and sample_id not in included:
            continue
        if sample_id in excluded:
            continue
        filtered.append(str(path))
    return filtered


# Identify ingroup FASTQs that live outside the project tree for transfer sync.
def external_ingroup_fastq_files():
    files = []
    for sample_id in INGROUP_STORAGE_SAMPLE_IDS:
        for read in ("1", "2"):
            fastq = ingroup_fastq_path(sample_id, read)
            if is_external_storage_path(fastq):
                files.append(fastq)
    return sorted(set(files))


outgroup_df = pd.read_csv(config["outgroups"], sep="\t")
OUTGROUP_SAMPLES = outgroup_df.groupby("sample_id")["srr_id"].apply(list).to_dict()
OUTGROUP_SAMPLE_IDS = list(OUTGROUP_SAMPLES.keys())
OUTGROUP_SRR_IDS = [srr for srrs in OUTGROUP_SAMPLES.values() for srr in srrs]

if "taxon" in outgroup_df.columns:
    OUTGROUP_TAXON_LABELS = (
        outgroup_df.assign(sample_id=outgroup_df["sample_id"].astype(str))
        .groupby("sample_id")["taxon"]
        .first()
        .apply(lambda value: re.sub(r"\s+", "_", str(value).strip()) if pd.notna(value) else None)
        .to_dict()
    )
else:
    OUTGROUP_TAXON_LABELS = {str(row.sample_id): None for row in outgroup_df.itertuples()}

OG_SEQ_COL = config.get("outgroups_seq_column", "sequencer")
LONGREAD_KEYWORDS = {keyword.lower() for keyword in config.get("longread_keywords", ["PacBio", "ONT", "Nanopore"])}


# Infer short-read vs long-read outgroups from the metadata sequencer field.
def _infer_read_type(sequencer):
    sequencer_name = str(sequencer or "").lower()
    for keyword in LONGREAD_KEYWORDS:
        if keyword in sequencer_name:
            return "long"
    return "short"


# Choose a minimap2 preset from the sequencer description.
def _minimap2_preset(sequencer):
    sequencer_name = str(sequencer or "").lower()
    if "ont" in sequencer_name or "nanopore" in sequencer_name:
        return "map-ont"
    if "hifi" in sequencer_name:
        return "map-hifi"
    return "map-pb"


if OG_SEQ_COL in outgroup_df.columns:
    OUTGROUP_READ_TYPE = (
        outgroup_df.groupby("sample_id")[OG_SEQ_COL].first().apply(_infer_read_type).to_dict()
    )
    OUTGROUP_MINIMAP2_PRESET = (
        outgroup_df.groupby("sample_id")[OG_SEQ_COL].first().apply(_minimap2_preset).to_dict()
    )
else:
    OUTGROUP_READ_TYPE = {sample_id: "short" for sample_id in OUTGROUP_SAMPLE_IDS}
    OUTGROUP_MINIMAP2_PRESET = {sample_id: "map-pb" for sample_id in OUTGROUP_SAMPLE_IDS}


# Resolve outgroup sample IDs from generic include/exclude/species selectors.
def outgroup_sample_ids(selected_ids=None, exclude_ids=None, species=None, include_all_when_unspecified=True):
    if selected_ids is None:
        sample_ids = list(OUTGROUP_SAMPLE_IDS) if include_all_when_unspecified else []
    else:
        sample_ids = [sample_id for sample_id in _parse_list(selected_ids) if sample_id in OUTGROUP_SAMPLE_IDS]

    species_list = _parse_species_list(species)
    if species_list:
        if "taxon" not in outgroup_df.columns:
            sample_ids = []
        else:
            allowed = set(
                outgroup_df.loc[
                    outgroup_df["taxon"].astype(str).isin(species_list),
                    "sample_id",
                ].astype(str)
            )
            sample_ids = [sample_id for sample_id in sample_ids if sample_id in allowed]

    excluded = set(_parse_list(exclude_ids))
    if excluded:
        sample_ids = [sample_id for sample_id in sample_ids if sample_id not in excluded]
    return sample_ids


# Analysis-friendly wrapper that returns no outgroups when the toggle is off.
def analysis_outgroup_sample_ids(include_outgroups, selected_ids=None, exclude_ids=None, species=None):
    if not include_outgroups:
        return []
    return outgroup_sample_ids(
        selected_ids=selected_ids,
        exclude_ids=exclude_ids,
        species=species,
        include_all_when_unspecified=True,
    )


# Materialize sliced outgroup BAM paths from selected sample IDs.
def outgroup_bam_paths(sample_ids):
    return [f"{OUTGROUP_SLICED_DIR}/{sample_id}.bam" for sample_id in sample_ids]


# Helper for Snakemake input sections that need sliced outgroup BAMs.
def sliced_outgroup_inputs(sample_ids):
    if not sample_ids:
        return []
    return expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=sample_ids)


# Write a plain one-BAM-per-line bamlist file.
def write_bamlist(output_path, bam_paths):
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(path) for path in bam_paths]
    output_file.write_text(("\n".join(lines) + "\n") if lines else "", encoding="ascii")


# Short taxon label used for TreeMix roots when taxon names are available.
def taxon_abbreviation(value):
    try:
        if pd.isna(value):
            return None
    except Exception:
        if value is None:
            return None
    words = str(value).strip().split()
    return "".join(word[:2] for word in words if word) or None


# Build a command-line outgroup option from the sample IDs present in a bamlist.
def outgroup_option_from_bamlist(bamlist_path, species=None, sample_ids=None):
    try:
        with open(bamlist_path) as handle:
            bamlist_ids = [bam_sample_id(line.strip()) for line in handle if line.strip()]
    except FileNotFoundError:
        return ""

    selected = sample_ids
    if selected is None:
        selected = outgroup_sample_ids(selected_ids=None, species=species, include_all_when_unspecified=True)
    selected_set = set(str(sample_id) for sample_id in selected)
    labels = [sample_id for sample_id in bamlist_ids if sample_id in selected_set]
    return "" if not labels else "--outgroup " + ",".join(labels)


# Build a TreeMix -root argument from selected outgroup samples.
def treemix_root_option(sample_ids, merge_outgroups=False, root_label=None):
    if not sample_ids:
        return ""
    og = outgroup_df[outgroup_df["sample_id"].astype(str).isin(sample_ids)].copy()
    labels = []
    if merge_outgroups:
        label = root_label
        if not label:
            if "taxon" in og.columns:
                abbreviations = [value for value in og["taxon"].apply(taxon_abbreviation).dropna().tolist() if value]
                unique = sorted(set(abbreviations))
                label = unique[0] if len(unique) == 1 else "OUT"
            else:
                label = "OUT"
        labels = [label]
    elif "taxon" in og.columns:
        labels = sorted(
            {
                value
                for value in og["taxon"].apply(taxon_abbreviation).dropna().tolist()
                if value
            }
        )
    else:
        labels = sorted(set(og["sample_id"].astype(str).tolist()))
    return "" if not labels else "-root " + ",".join(labels)
