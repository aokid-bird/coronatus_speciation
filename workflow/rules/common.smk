import pandas as pd
from pathlib import Path
from itertools import combinations
import re
import os
from datetime import datetime, timezone

_CFG_MISSING = object()


def _config_value(mapping, *keys, default=_CFG_MISSING):
    value = mapping
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def _set_config_value(mapping, keys, value):
    target = mapping
    for key in keys[:-1]:
        child = target.get(key)
        if not isinstance(child, dict):
            child = {}
            target[key] = child
        target = child
    target[keys[-1]] = value


def _backfill_config(target_keys, *source_options, default=_CFG_MISSING):
    current = _config_value(config, *target_keys, default=_CFG_MISSING)
    if current is not _CFG_MISSING:
        return current
    for source_keys in source_options:
        candidate = _config_value(config, *source_keys, default=_CFG_MISSING)
        if candidate is not _CFG_MISSING:
            _set_config_value(config, target_keys, candidate)
            return candidate
    if default is not _CFG_MISSING:
        _set_config_value(config, target_keys, default)
        return default
    return _CFG_MISSING


# Backfill legacy top-level config keys from the grouped refactor schema so
# existing rule files can remain stable while defaults move to clearer sections.
_backfill_config(("environment",), ("project", "environment"), default="cluster")
_backfill_config(("threads",), ("project", "threads"), default=1)
_backfill_config(("output_prefix",), ("project", "output_prefix"), default="defaults")
_backfill_config(("populations",), ("project", "populations"), default=[])
_backfill_config(("group_col",), ("project", "group_col"), default=None)
_backfill_config(("population_labels",), ("project", "population_labels"), default=None)

_backfill_config(("bam_dir",), ("paths", "bam_dir"), default="results/bwa")
_backfill_config(("singularity_dir",), ("paths", "singularity_dir"), default="${HOME}/envs/singularity")
_backfill_config(("transfer",), ("paths", "transfer"), default={})

_backfill_config(("samples",), ("inputs", "samples_tsv"), default="data/samples.tsv")
_backfill_config(("outgroups",), ("inputs", "outgroups_tsv"), default="data/outgroup.tsv")
_backfill_config(("references_tsv",), ("inputs", "references_tsv"), default="data/references.tsv")

_backfill_config(("reference_download_method",), ("reference", "download_method"), default="datasets")
_backfill_config(("reference_dir",), ("reference", "storage_dir"), default="data/reference")
_backfill_config(("outgroups_seq_column",), ("metadata", "outgroups_seq_column"), default="sequencer")
_backfill_config(("longread_keywords",), ("metadata", "longread_keywords"), default=["PacBio", "ONT", "Nanopore"])

_backfill_config(("angsd_common_args",), ("angsd", "common_args"), default="")
_backfill_config(("angsd_args",), ("angsd", "args"), default={})
_backfill_config(("minIndRatio",), ("angsd", "min_ind_ratio"), default={})

for analysis_name in (
    "angsd_intersect",
    "angsd_global",
    "angsd_global_unrelated_unlinked",
    "angsd_sfs",
    "angsd_raxml",
    "slice_outgroups",
    "ngsrelate",
    "ngsld",
    "pcangsd",
    "ngsadmix",
    "ngsdist",
    "raxml",
    "treemix",
    "abbababa2",
    "snapp",
):
    _backfill_config((analysis_name,), ("analyses", analysis_name), default={})

_backfill_config(("sfs_analysis",), ("analyses", "sfs"), default={})

# main parameters
output_prefix = config["output_prefix"] # scenario symbols
groups = config["populations"]
group_col = config["group_col"]
LEGACY_GLOBAL_THREADS = int(config.get("threads", 1))


def _config_section(*keys):
    value = config
    for key in keys:
        if not isinstance(value, dict):
            return {}
        value = value.get(key, {})
    return value if isinstance(value, dict) else {}


def _resolve_threads(section, default, legacy_fallback=True):
    value = section.get("threads")
    if value is None and legacy_fallback:
        value = config.get("threads", default)
    if value is None:
        value = default
    return int(value)

_POP_LABEL_CFG = config.get("population_labels") or {}
if isinstance(_POP_LABEL_CFG, dict):
    POPULATION_LABELS = {str(k): str(v) for k, v in _POP_LABEL_CFG.items()}
else:
    POPULATION_LABELS = {}

def get_population_label(pop: str) -> str:
    return POPULATION_LABELS.get(str(pop), str(pop))


def _is_enabled(section, default=True):
    return bool(section.get("enabled", default))

SFS_CFG = (config.get("sfs_analysis", {}) or {})

def _norm_token(value, default):
    v = str(value).strip().lower().replace(" ", "_") if value is not None else None
    return v if v else default

SFS_SITE_FILTERS = [_norm_token(v, "linked") for v in SFS_CFG.get("site_filters", ["linked"]) or ["linked"]]
SFS_FOLD_STATES = [_norm_token(v, "fold") for v in SFS_CFG.get("fold_states", ["fold"]) or ["fold"]]
SFS_INCLUDE_MONOMORPHIC = bool(SFS_CFG.get("include_monomorphic", True))
_SFS_REAL_CFG = (SFS_CFG.get("realSFS", {}) or {})
SFS_MAXITER = int(_SFS_REAL_CFG.get("maxiter", 50000))
SFS_TOLE = str(_SFS_REAL_CFG.get("tole", "1e-6"))
SFS_NEEDS_UNLINKED_SITES = "unlinked" in {v.lower() for v in SFS_SITE_FILTERS}
SFS_PAIRWISE_COMBOS = list(combinations(groups, 2)) if len(groups) >= 2 else []

# directories
STORAGE_CFG = config.get("storage", {}) or {}

def _storage_path(*keys, default=None):
    value = STORAGE_CFG
    for key in keys:
        if not isinstance(value, dict):
            return default
        value = value.get(key)
    return default if value in (None, "") else value

config_bam_base = config["bam_dir"]
config_bam_dir = _storage_path("bam", "ingroup_dir", default=f"{config_bam_base}/{output_prefix}")
OUTGROUP_BAM_DIR = _storage_path("bam", "outgroup_dir", default=f"{config_bam_dir}/outgroups")
OUTGROUP_SLICED_DIR = f"results/outgroups_sliced/{output_prefix}"
TRANSFER_CFG = config.get("transfer", {}) or {}
CLUSTER_STORAGE_ROOT_DEFAULT = str(TRANSFER_CFG.get("cluster_storage_root", "/lfs/aokid"))
PROJECT_ROOT = Path.cwd().resolve()
config_singularity_dir = os.path.expandvars(config["singularity_dir"])
# other parameters
kin_thr = config["ngsrelate"]["kinship_threshold"]

#=================#
#### REFERENCE ####
#=================#
def derive_indexes(fasta: str):
    p = Path(fasta)
    # adaptation to fasta file variations (e.g.,.fa.gz) 
    if p.suffix == ".gz":
        # Keep the .gz extension in the basename so samtools faidx outputs match
        fai = f"{p}.fai"
        gzi = f"{p}.gzi"
    else:
        fai = f"{p}.fai"
        gzi = None
    # # .dict（(if needed))
    # dict_path = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", ".dict", str(p))
    # return fai, gzi, dict_path
    return fai, gzi

ENVIRONMENT = config.get("environment", "cluster")

# from config
REF_CONFIG_PATH = config.get("reference", {}).get("fasta", config.get("ref"))
REFERENCES_TSV = config.get("references_tsv", "data/references.tsv")
REFERENCE_DIR = _storage_path("reference", "dir", default=config.get("reference_dir", "data/reference"))
REFERENCE_METHOD = config.get("reference_download_method", "datasets")
OUTGROUP_RAW_DIR = _storage_path("outgroup", "raw_dir", default="data/raw/outgroup")
OUTGROUP_MERGED_DIR = _storage_path("outgroup", "merged_dir", default="data/merged/outgroup")
INGROUP_TRIM_DIR = _storage_path("derived", "ingroup", "trim_dir", default="results/trimmomatic/ingroup")
OUTGROUP_TRIM_DIR = _storage_path("derived", "outgroup", "trim_dir", default="results/trimmomatic/outgroup")
INGROUP_QC_BASE_DIR = _storage_path("derived", "ingroup", "qc_dir", default="results/qc/ingroup")
OUTGROUP_QC_BASE_DIR = _storage_path("derived", "outgroup", "qc_dir", default="results/qc/outgroup")
INGROUP_MAP_TMP_DIR = _storage_path("derived", "ingroup", "mapping_tmp_dir", default=f"results/mapping/{output_prefix}/ingroup")
OUTGROUP_MAP_TMP_DIR = _storage_path("derived", "outgroup", "mapping_tmp_dir", default=f"results/mapping/{output_prefix}/outgroup")
OUTGROUP_LR_FILTER_DIR = _storage_path("derived", "outgroup", "longread_filter_dir", default="results/longread/filter")

def manifest_paths(directory: str):
    return (
        str(Path(directory) / "README.md"),
        str(Path(directory) / "provenance.yaml"),
    )

def _is_directory_target_writable(directory: str) -> bool:
    probe = Path(directory)
    while not probe.exists():
        if probe.parent == probe:
            return False
        probe = probe.parent
    return probe.is_dir() and os.access(probe, os.W_OK)

def manifest_directory(directory: str, fallback_directory: str = None) -> str:
    if _is_directory_target_writable(directory):
        return directory
    if fallback_directory:
        return fallback_directory
    return directory

REFERENCE_MANIFEST_DIR = manifest_directory(
    REFERENCE_DIR,
    fallback_directory=f"results/reference/{output_prefix}/source_metadata",
)

def write_storage_manifest(directory: str, title: str, producer: str, details=None):
    details = details or {}
    directory_path = Path(directory)
    directory_path.mkdir(parents=True, exist_ok=True)
    readme_path = directory_path / "README.md"
    yaml_path = directory_path / "provenance.yaml"
    timestamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    project_dir = str(Path.cwd())
    pipeline_name = Path(project_dir).name

    yaml_lines = [
        f'title: "{title}"',
        f'directory: "{directory_path}"',
        f'created_utc: "{timestamp}"',
        f'producer: "{producer}"',
        f'pipeline_name: "{pipeline_name}"',
        f'project_dir: "{project_dir}"',
        f'output_prefix: "{output_prefix}"',
    ]
    if details:
        yaml_lines.append("details:")
        for key, value in details.items():
            safe_key = str(key).replace(":", "_")
            safe_value = str(value).replace('"', '\\"')
            yaml_lines.append(f'  {safe_key}: "{safe_value}"')

    readme_lines = [
        f"# {title}",
        "",
        "This directory is managed by the `gbs_pipeline` Snakemake workflow.",
        "",
        f"- Directory: `{directory_path}`",
        f"- Created (UTC): `{timestamp}`",
        f"- Producer: `{producer}`",
        f"- Pipeline: `{pipeline_name}`",
        f"- Project: `{project_dir}`",
        f"- output_prefix at write time: `{output_prefix}`",
    ]
    if details:
        readme_lines.extend(["", "## Details"])
        for key, value in details.items():
            readme_lines.append(f"- {key}: `{value}`")

    readme_path.write_text("\n".join(readme_lines) + "\n", encoding="ascii")
    yaml_path.write_text("\n".join(yaml_lines) + "\n", encoding="ascii")

# metadata-driven reference (only on local)
REF_FROM_METADATA = False
REFERENCE_ACCESSION = None
REFERENCE_NAME = None
REFERENCE_URL = ""
REFERENCE_ZIP = None

if REF_CONFIG_PATH:
    REF = REF_CONFIG_PATH
else:
    if ENVIRONMENT == "local":
        # read references.tsv (columns: accession, name, optional url)
        ref_df = pd.read_csv(REFERENCES_TSV, sep="\t")
        if ref_df.empty:
            raise ValueError(f"references.tsv has no rows: {REFERENCES_TSV}")
        row = ref_df.iloc[0]
        REFERENCE_ACCESSION = str(row.get("accession"))
        REFERENCE_NAME = str(row.get("name", REFERENCE_ACCESSION))
        REFERENCE_URL = str(row.get("url", ""))
        Path(f"{REFERENCE_DIR}/raw").mkdir(parents=True, exist_ok=True)
        REFERENCE_ZIP = f"{REFERENCE_DIR}/raw/{REFERENCE_ACCESSION}.zip"
        REF = f"{REFERENCE_DIR}/{REFERENCE_NAME}.fna.gz"
        REF_FROM_METADATA = True
    else:
        # On cluster, user is expected to provide reference path in config
        raise ValueError(
            "On cluster environment, set reference.fasta in config to the transferred reference."
        )

# index paths and auxiliary files
FAI_DEFAULT, GZI = derive_indexes(REF)
FAI  = config.get("reference", {}).get("fai") or FAI_DEFAULT
REF_CHR = f"{REFERENCE_DIR}/{Path(REF).stem}.chr" if REF_FROM_METADATA else f"{Path(REF).with_suffix('').with_suffix('').as_posix()}.chr"

#===============#
#### MAPPING ####
#===============#
# Ingroup sample metadata
_SAMPLES_DF = pd.read_csv(config["samples"], sep="\t")
INGROUP_SAMPLE_IDS = _SAMPLES_DF["sample"].astype(str).tolist()
READS_CFG = (config.get("reads", {}) or {})
INGROUP_READS_DIR = READS_CFG.get("ingroup_dir", "data/raw/ingroup")
INGROUP_READS_META_CFG = (READS_CFG.get("ingroup_metadata", {}) or {})

INGROUP_FASTQ_DIR_COL = str(INGROUP_READS_META_CFG.get("dir_col", "fastq_dir"))
INGROUP_FASTQ_PREFIX_COL = str(INGROUP_READS_META_CFG.get("prefix_col", "fastq_prefix"))
INGROUP_FASTQ_R1_SUFFIX_COL = str(INGROUP_READS_META_CFG.get("r1_suffix_col", "fastq_r1_suffix"))
INGROUP_FASTQ_R2_SUFFIX_COL = str(INGROUP_READS_META_CFG.get("r2_suffix_col", "fastq_r2_suffix"))
INGROUP_FASTQ_EXT_COL = str(INGROUP_READS_META_CFG.get("extension_col", "fastq_extension"))

INGROUP_FASTQ_DEFAULT_PREFIX = str(INGROUP_READS_META_CFG.get("default_prefix", ""))
INGROUP_FASTQ_DEFAULT_R1_SUFFIX = str(INGROUP_READS_META_CFG.get("default_r1_suffix", "_1"))
INGROUP_FASTQ_DEFAULT_R2_SUFFIX = str(INGROUP_READS_META_CFG.get("default_r2_suffix", "_2"))
INGROUP_FASTQ_DEFAULT_EXT = str(INGROUP_READS_META_CFG.get("default_extension", ".fastq.gz"))

def _sample_value(row, colname, default=""):
    if colname in row.index and pd.notna(row[colname]):
        return str(row[colname])
    return default

def _normalize_fastq_ext(ext: str) -> str:
    value = str(ext or "").strip()
    if not value:
        return ""
    return value if value.startswith(".") else f".{value}"

_INGROUP_FASTQ_RECORDS = {}
for _, _row in _SAMPLES_DF.iterrows():
    sample_id = str(_row["sample"])
    sample_dir = _sample_value(_row, INGROUP_FASTQ_DIR_COL, INGROUP_READS_DIR).strip()
    prefix = _sample_value(_row, INGROUP_FASTQ_PREFIX_COL, INGROUP_FASTQ_DEFAULT_PREFIX)
    r1_suffix = _sample_value(_row, INGROUP_FASTQ_R1_SUFFIX_COL, INGROUP_FASTQ_DEFAULT_R1_SUFFIX)
    r2_suffix = _sample_value(_row, INGROUP_FASTQ_R2_SUFFIX_COL, INGROUP_FASTQ_DEFAULT_R2_SUFFIX)
    ext = _normalize_fastq_ext(_sample_value(_row, INGROUP_FASTQ_EXT_COL, INGROUP_FASTQ_DEFAULT_EXT))
    if not sample_dir:
        raise ValueError(f"Sample '{sample_id}' has an empty FASTQ directory after config/metadata resolution.")
    _INGROUP_FASTQ_RECORDS[sample_id] = {
        "dir": sample_dir,
        "prefix": prefix,
        "r1_suffix": r1_suffix,
        "r2_suffix": r2_suffix,
        "extension": ext,
    }

def ingroup_fastq_path(sample_id: str, read: str) -> str:
    sid = str(sample_id)
    if sid not in _INGROUP_FASTQ_RECORDS:
        raise KeyError(f"Unknown ingroup sample_id '{sid}'")
    if str(read) not in {"1", "2"}:
        raise ValueError(f"read must be '1' or '2', got '{read}'")
    record = _INGROUP_FASTQ_RECORDS[sid]
    suffix = record["r1_suffix"] if str(read) == "1" else record["r2_suffix"]
    filename = f"{record['prefix']}{sid}{suffix}{record['extension']}"
    return str(Path(record["dir"]) / filename)


def bam_sample_id(path: str) -> str:
    name = Path(str(path)).name
    for suffix in (".bam", ".cram", ".sam"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if "_slice" in name:
        name = name.split("_slice", 1)[0]
    return name


def ingroup_sample_ids(populations=None, exclude_samples=None):
    df = _SAMPLES_DF
    if populations:
        df = df[df[group_col].isin(populations)]
    exclude = set(_parse_exclude_list(exclude_samples))
    sample_ids = df["sample"].astype(str).tolist()
    if exclude:
        sample_ids = [sid for sid in sample_ids if sid not in exclude]
    return sample_ids


def ingroup_bam_paths(sample_ids=None, populations=None, exclude_samples=None, bam_dir=None):
    selected = list(sample_ids) if sample_ids is not None else ingroup_sample_ids(
        populations=populations,
        exclude_samples=exclude_samples,
    )
    target_dir = bam_dir or config_bam_dir
    return [f"{target_dir}/{sid}.bam" for sid in selected]


def filter_bam_paths_by_sample_ids(bam_paths, exclude_samples=None, include_samples=None):
    excluded = set(_parse_exclude_list(exclude_samples))
    included = None if include_samples is None else set(str(s) for s in include_samples)
    filtered = []
    for path in bam_paths:
        sample_id = bam_sample_id(path)
        if included is not None and sample_id not in included:
            continue
        if sample_id in excluded:
            continue
        filtered.append(str(path))
    return filtered

def _path_is_within(path: Path, base: Path) -> bool:
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False

def is_external_storage_path(path: str) -> bool:
    candidate = Path(str(path)).expanduser()
    if not candidate.is_absolute():
        return False
    return not _path_is_within(candidate, PROJECT_ROOT)

def remote_mirror_path(path: str, cluster_root: str = None) -> str:
    root = Path(cluster_root or CLUSTER_STORAGE_ROOT_DEFAULT)
    candidate = Path(str(path)).expanduser()
    return str(root / str(candidate).lstrip("/"))

def external_ingroup_fastq_files():
    files = []
    for sample_id in INGROUP_SAMPLE_IDS:
        for read in ("1", "2"):
            fq = ingroup_fastq_path(sample_id, read)
            if is_external_storage_path(fq):
                files.append(fq)
    return sorted(set(files))

def external_transfer_directories():
    candidates = [
        REFERENCE_DIR,
        OUTGROUP_RAW_DIR,
        OUTGROUP_MERGED_DIR,
        INGROUP_TRIM_DIR,
        OUTGROUP_TRIM_DIR,
        INGROUP_QC_BASE_DIR,
        OUTGROUP_QC_BASE_DIR,
        OUTGROUP_LR_FILTER_DIR,
        config_bam_dir,
        OUTGROUP_BAM_DIR,
    ]
    return sorted(set(path for path in candidates if is_external_storage_path(path)))

def _merge_mapping_scope(scope: str):
    base = {k: v for k, v in (config.get("mapping", {}) or {}).items() if k not in {"ingroup", "outgroup"}}
    scoped = (config.get("mapping", {}) or {}).get(scope, {})
    if isinstance(scoped, dict):
        merged = base.copy()
        merged.update(scoped)
        return merged
    return base

MAPPING_INGROUP_CFG = _merge_mapping_scope("ingroup")
MAPPING_OUTGROUP_CFG = _merge_mapping_scope("outgroup")
MAPPING_INGROUP_THREADS = _resolve_threads(MAPPING_INGROUP_CFG, 6)
MAPPING_OUTGROUP_THREADS = _resolve_threads(MAPPING_OUTGROUP_CFG, 6)

# Choose mapping backend per scope: bwa or bwa-mem2
MAPPER_INGROUP = MAPPING_INGROUP_CFG.get("mapper", config.get("mapper", "bwa"))
MAPPER_OUTGROUP = MAPPING_OUTGROUP_CFG.get("mapper", MAPPER_INGROUP)

# Uncompressed FASTA path for mappers that require it
ref_path = Path(REF)
if ref_path.suffix == ".gz":
    REF_UNZIPPED = str(ref_path.with_suffix(""))
else:
    REF_UNZIPPED = str(ref_path)

# Optional separate directory to store mapper indices (useful if reference dir is read-only)
REF_INDEX_DIR = (config.get("reference", {}) or {}).get("index_dir", None)
if REF_INDEX_DIR:
    REF_INDEX_PREFIX_BASE = str(Path(REF_INDEX_DIR) / Path(REF_UNZIPPED).name)
else:
    REF_INDEX_PREFIX_BASE = REF_UNZIPPED

# Use distinct index prefixes when bwa and bwa-mem2 are both requested to avoid clobbering files.
REF_INDEX_PREFIX_BWA = REF_INDEX_PREFIX_BASE
REF_INDEX_PREFIX_BWAMEM2 = REF_INDEX_PREFIX_BASE
if MAPPER_INGROUP != MAPPER_OUTGROUP:
    REF_INDEX_PREFIX_BWA = f"{REF_INDEX_PREFIX_BASE}.bwa"
    REF_INDEX_PREFIX_BWAMEM2 = f"{REF_INDEX_PREFIX_BASE}.bwamem2"

# Representative index artifact per mapper (used as rule output) and mapper arg
BWA_INDEX_MAIN = f"{REF_INDEX_PREFIX_BWA}.bwt"
BWAMEM2_INDEX_MAIN = f"{REF_INDEX_PREFIX_BWAMEM2}.0123"
# Index selection per mapper scope
def _ref_index_for(mapper: str):
    return BWAMEM2_INDEX_MAIN if mapper == "bwa-mem2" else BWA_INDEX_MAIN

REF_MAP_INDEX_INGROUP = _ref_index_for(MAPPER_INGROUP)
REF_MAP_INDEX_OUTGROUP = _ref_index_for(MAPPER_OUTGROUP)
REF_MAP_INDEX_MAIN = REF_MAP_INDEX_INGROUP
REF_MAP_INDEXES = sorted({REF_MAP_INDEX_INGROUP, REF_MAP_INDEX_OUTGROUP})

# Argument to pass to mapper (prefix if indexed elsewhere; equals REF_UNZIPPED if next to fasta)
REF_MAP_ARG_BWA = REF_INDEX_PREFIX_BWA
REF_MAP_ARG_BWAMEM2 = REF_INDEX_PREFIX_BWAMEM2
REF_MAP_ARG_INGROUP = REF_MAP_ARG_BWAMEM2 if MAPPER_INGROUP == "bwa-mem2" else REF_MAP_ARG_BWA
REF_MAP_ARG_OUTGROUP = REF_MAP_ARG_BWAMEM2 if MAPPER_OUTGROUP == "bwa-mem2" else REF_MAP_ARG_BWA
REF_MAP_ARG = REF_MAP_ARG_INGROUP
REF_INDEX_PREFIX = REF_INDEX_PREFIX_BASE  # backwards compatibility
MAPPER = MAPPER_INGROUP  # backwards compatibility


def _parse_exclude_list(raw):
    if raw is None:
        return []
    if isinstance(raw, str):
        return [s.strip() for s in raw.split(",") if s.strip()]
    return [str(s).strip() for s in raw if str(s).strip()]


def _parse_species_list(raw):
    return _parse_exclude_list(raw)

#================#
#### TREEMIX  ####
#================#
TREEMIX_CFG = config.get("treemix", {}) or {}
GLACTOOLS_BIN = TREEMIX_CFG.get("glactools_bin", "workflow/bin/glactools")
TREEMIX_MODE = TREEMIX_CFG.get("mode", "pop")  # pop or indiv
TREEMIX_INCLUDE_OUTGROUPS = bool(TREEMIX_CFG.get("include_outgroups", True))
TREEMIX_ROOT_LABEL = TREEMIX_CFG.get("root_label")
TREEMIX_MAX_M = int(TREEMIX_CFG.get("max_m", 6))
TREEMIX_REPS = int(TREEMIX_CFG.get("reps", 10))
TREEMIX_THREADS = _resolve_threads(TREEMIX_CFG, LEGACY_GLOBAL_THREADS)
TREEMIX_TIMEOUT_SECONDS = int(TREEMIX_CFG.get("timeout_seconds", 300))
TREEMIX_MAX_ATTEMPTS = int(TREEMIX_CFG.get("max_attempts", 3))
# Optional list of sample IDs to exclude from TreeMix analysis
_treemix_exclude_cfg = TREEMIX_CFG.get("exclude_samples") or []
if isinstance(_treemix_exclude_cfg, str):
    TREEMIX_EXCLUDE_SAMPLES = [s.strip() for s in _treemix_exclude_cfg.split(",") if s.strip()]
else:
    TREEMIX_EXCLUDE_SAMPLES = [str(s).strip() for s in _treemix_exclude_cfg if str(s).strip()]
_treemix_exclude_out_cfg = TREEMIX_CFG.get("exclude_outgroups") or []
if isinstance(_treemix_exclude_out_cfg, str):
    TREEMIX_EXCLUDE_OUTGROUPS = [s.strip() for s in _treemix_exclude_out_cfg.split(",") if s.strip()]
else:
    TREEMIX_EXCLUDE_OUTGROUPS = [str(s).strip() for s in _treemix_exclude_out_cfg if str(s).strip()]
TREEMIX_OUTGROUP_SPECIES = _parse_species_list(TREEMIX_CFG.get("outgroup_species"))
# plotting_funcs: prefer fixed repo path if present; no config knob
_PLOTFUNC_DEFAULT = Path("workflow/scripts/treemix_plotting_funcs.R")
TREEMIX_PLOTTING_FUNCS = str(_PLOTFUNC_DEFAULT) if _PLOTFUNC_DEFAULT.exists() else None

#===========================#
#### ANALYSIS OUTGROUPS  ####
#===========================#
# Per-analysis toggles to include outgroups (sliced BAMs)
ANGSD_GLOBAL_CFG = (config.get("angsd_global", {}) or {})
ANGSD_RAXML_CFG  = (config.get("angsd_raxml", {}) or {})
ANGSD_INTERSECT_CFG = (config.get("angsd_intersect", {}) or {})
ANGSD_SFS_CFG = (config.get("angsd_sfs", {}) or {})
ANGSD_SNAPP_CFG = (config.get("snapp", {}) or {})
ANGSD_GLOBAL_INCLUDE_OUTGROUPS = bool(ANGSD_GLOBAL_CFG.get("include_outgroups", False))
ANGSD_RAXML_INCLUDE_OUTGROUPS  = bool(ANGSD_RAXML_CFG.get("include_outgroups", False))
NGSDIST_CFG = _config_section("ngsdist")
NGSDIST_INCLUDE_OUTGROUPS = bool(NGSDIST_CFG.get("include_outgroups", ANGSD_GLOBAL_INCLUDE_OUTGROUPS))
ANGSD_INTERSECT_THREADS = _resolve_threads(ANGSD_INTERSECT_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_GLOBAL_THREADS = _resolve_threads(ANGSD_GLOBAL_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_GLOBAL_UNRELATED_THREADS = _resolve_threads(
    _config_section("angsd_global_unrelated_unlinked"),
    ANGSD_GLOBAL_THREADS,
)
ANGSD_SFS_THREADS = _resolve_threads(ANGSD_SFS_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_RAXML_THREADS = _resolve_threads(ANGSD_RAXML_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_SNAPP_THREADS = _resolve_threads(ANGSD_SNAPP_CFG, ANGSD_GLOBAL_THREADS)
SLICE_OUTGROUPS_THREADS = _resolve_threads(_config_section("slice_outgroups"), 6)
NGSLD_THREADS = _resolve_threads(_config_section("ngsld"), LEGACY_GLOBAL_THREADS)
REALSFS_THREADS = _resolve_threads(_config_section("sfs_analysis", "realSFS"), 10, legacy_fallback=False)
QC_INGROUP_THREADS = _resolve_threads(_config_section("qc", "ingroup", "trimmomatic"), 4)
QC_OUTGROUP_THREADS = _resolve_threads(_config_section("qc", "outgroup", "trimmomatic"), 4)

ANGSD_RAXML_DOWNSAMPLE_CFG = (ANGSD_RAXML_CFG.get("downsampling") or {})

def _parse_max_per_population(raw):
    if raw is None or (isinstance(raw, str) and raw.strip().lower() in {"", "none"}):
        return None
    if isinstance(raw, dict):
        parsed = {}
        for key, value in raw.items():
            if value is None or (isinstance(value, str) and value.strip().lower() in {"", "none"}):
                parsed[str(key)] = None
            else:
                parsed[str(key)] = int(value)
        return parsed
    return int(raw)

def _parse_optional_list(raw):
    if raw is None:
        return None
    parsed = _parse_exclude_list(raw)
    return parsed if parsed else []

def _normalise_seed(value):
    if value is None or (isinstance(value, str) and value.strip().lower() in {"", "none"}):
        return None
    return int(value)

ANGSD_RAXML_DOWNSAMPLE_MAX = None
try:
    ANGSD_RAXML_DOWNSAMPLE_MAX = _parse_max_per_population(
        ANGSD_RAXML_DOWNSAMPLE_CFG.get("max_per_population")
    )
except (TypeError, ValueError):
    ANGSD_RAXML_DOWNSAMPLE_MAX = None

ANGSD_RAXML_DOWNSAMPLE_EXCLUDE = _parse_exclude_list(
    ANGSD_RAXML_DOWNSAMPLE_CFG.get("exclude_samples")
)
ANGSD_RAXML_DOWNSAMPLE_USE_ALL = bool(ANGSD_RAXML_DOWNSAMPLE_CFG.get("use_all_samples", False))
_seed_raw = ANGSD_RAXML_DOWNSAMPLE_CFG.get("seed", None)
try:
    ANGSD_RAXML_DOWNSAMPLE_SEED = _normalise_seed(_seed_raw)
except (TypeError, ValueError):
    ANGSD_RAXML_DOWNSAMPLE_SEED = None

# Helper: resolve minIndRatio for a given key with fallbacks
def get_minInd_ratio(key, default=None):
    m = (config.get("minIndRatio", {}) or {})
    if key in m:
        return float(m[key])
    # common fallback to 'global' or 'intersect'
    if key != "global" and "global" in m:
        return float(m["global"])
    if key != "intersect" and "intersect" in m:
        return float(m["intersect"])
    return default

#================#
#### OUTGROUP ####
#================#
# Read outgroup SRR info from TSV
outgroup_df = pd.read_csv(config["outgroups"], sep="\t")

# Group SRR IDs by BioSample names
OUTGROUP_SAMPLES = (
    outgroup_df
    .groupby("sample_id")["srr_id"]
    .apply(list)
    .to_dict()
)
OUTGROUP_SAMPLE_IDS = list(OUTGROUP_SAMPLES.keys())

if "taxon" in outgroup_df.columns:
    OUTGROUP_TAXON_LABELS = (
        outgroup_df.assign(sample_id=outgroup_df["sample_id"].astype(str))
        .groupby("sample_id")["taxon"]
        .first()
        .apply(lambda x: re.sub(r"\s+", "_", str(x).strip()) if pd.notna(x) else None)
        .to_dict()
    )
else:
    OUTGROUP_TAXON_LABELS = {str(row.sample_id): None for row in outgroup_df.itertuples()}

# Flattened list of all SRR IDs across samples (for fetch_sra targets)
OUTGROUP_SRR_IDS = [srr for srrs in OUTGROUP_SAMPLES.values() for srr in srrs]


def outgroup_sample_ids(selected_ids=None, exclude_ids=None, species=None, include_all_when_unspecified=True):
    if selected_ids is None:
        sample_ids = list(OUTGROUP_SAMPLE_IDS) if include_all_when_unspecified else []
    else:
        sample_ids = [sid for sid in _parse_exclude_list(selected_ids) if sid in OUTGROUP_SAMPLE_IDS]

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
            sample_ids = [sid for sid in sample_ids if sid in allowed]

    exclude = set(_parse_exclude_list(exclude_ids))
    if exclude:
        sample_ids = [sid for sid in sample_ids if sid not in exclude]
    return sample_ids


def analysis_outgroup_sample_ids(include_outgroups, selected_ids=None, exclude_ids=None, species=None):
    if not include_outgroups:
        return []
    return outgroup_sample_ids(
        selected_ids=selected_ids,
        exclude_ids=exclude_ids,
        species=species,
        include_all_when_unspecified=True,
    )


def outgroup_bam_paths(sample_ids):
    return [f"{OUTGROUP_SLICED_DIR}/{sid}.bam" for sid in sample_ids]


def sliced_outgroup_inputs(sample_ids):
    if not sample_ids:
        return []
    return expand(f"{OUTGROUP_SLICED_DIR}/{{sample_id}}.bam", sample_id=sample_ids)


def write_bamlist(output_path: str, bam_paths):
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(path) for path in bam_paths]
    output_file.write_text(
        ("\n".join(lines) + "\n") if lines else "",
        encoding="ascii",
    )


# Determine read type per outgroup sample (short vs long) from sequencer column
OG_SEQ_COL = config.get("outgroups_seq_column", "sequencer")
LONGREAD_KEYWORDS = set(config.get("longread_keywords", ["PacBio", "ONT", "Nanopore"]))

def _infer_read_type(sequencer: str) -> str:
    s = str(sequencer or "").lower()
    for kw in LONGREAD_KEYWORDS:
        if kw.lower() in s:
            return "long"
    return "short"

if OG_SEQ_COL in outgroup_df.columns:
    OUTGROUP_READ_TYPE = (
        outgroup_df.groupby("sample_id")[OG_SEQ_COL].first().apply(_infer_read_type).to_dict()
    )
else:
    # default all to short if no column provided
    OUTGROUP_READ_TYPE = {sid: "short" for sid in OUTGROUP_SAMPLE_IDS}

# Minimap2 preset for long-read samples (simple heuristic). Can be overridden per sample via config.
def _minimap2_preset(sequencer: str) -> str:
    s = str(sequencer or "").lower()
    if "ont" in s or "nanopore" in s:
        return "map-ont"
    if "hifi" in s:
        return "map-hifi"
    # default PacBio CLR/unknown long reads
    return "map-pb"

if OG_SEQ_COL in outgroup_df.columns:
    OUTGROUP_MINIMAP2_PRESET = (
        outgroup_df.groupby("sample_id")[OG_SEQ_COL].first().apply(_minimap2_preset).to_dict()
    )
else:
    OUTGROUP_MINIMAP2_PRESET = {sid: "map-pb" for sid in OUTGROUP_SAMPLE_IDS}


ABBABABA2_CFG = (config.get("abbababa2", {}) or {})
ABBABABA2_ENABLED = bool(ABBABABA2_CFG.get("enabled", False))
ABBABABA2_THREADS = int(ABBABABA2_CFG.get("threads", LEGACY_GLOBAL_THREADS))
_abb_args = str(ABBABABA2_CFG.get("angsd_args", "")).strip()
if not _abb_args:
    _abb_args = "-doAbbababa2 1 -doCounts 1 -minMapQ 30 -minQ 20 -baq 2 -useLast 1"
ABBABABA2_ANGSD_ARGS = _abb_args

ABBABABA2_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=True,
    selected_ids=ABBABABA2_CFG.get("outgroup_samples"),
)

_abb_label_cfg = ABBABABA2_CFG.get("outgroup_label")
ABBABABA2_OUTGROUP_LABEL_DEFAULT = None
if isinstance(_abb_label_cfg, dict):
    ABBABABA2_OUTGROUP_LABELS = {
        str(k): str(v)
        for k, v in _abb_label_cfg.items()
        if str(v).strip()
    }
else:
    ABBABABA2_OUTGROUP_LABELS = {}
    if _abb_label_cfg is not None:
        _lbl = str(_abb_label_cfg).strip()
        if _lbl:
            ABBABABA2_OUTGROUP_LABEL_DEFAULT = _lbl

_abb_exclude = ABBABABA2_CFG.get("exclude_samples") or []
if isinstance(_abb_exclude, str):
    ABBABABA2_EXCLUDE_SAMPLES = [s.strip() for s in _abb_exclude.split(",") if s.strip()]
else:
    ABBABABA2_EXCLUDE_SAMPLES = [str(s).strip() for s in _abb_exclude if str(s).strip()]

def resolve_abbababa2_outgroup_label(sample_id: str) -> str:
    sid = str(sample_id)
    if sid in ABBABABA2_OUTGROUP_LABELS:
        return ABBABABA2_OUTGROUP_LABELS[sid]
    if ABBABABA2_OUTGROUP_LABEL_DEFAULT:
        return ABBABABA2_OUTGROUP_LABEL_DEFAULT
    fallback = OUTGROUP_TAXON_LABELS.get(sid)
    if fallback:
        return fallback
    return sid


rule bcf2vcf:
    """
    Convert bcf to vcf using bcftools
    """
    input:
        bcf = lambda wc: f"results/angsd_{wc.angsd_run}/{output_prefix}/gl.bcf"
    output:
        vcf = f"results/angsd_{{angsd_run}}/{output_prefix}/gl.vcf.gz"
    conda:
        "../envs/bcftools_env.yaml"
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf}
        """

#================#
####  SNAPP   ####
#================#
SNAPP_CFG = (config.get("snapp", {}) or {})
SNAPP_ENABLED = bool(SNAPP_CFG.get("enabled", False))
SNAPP_MAX_PER_POP = int(SNAPP_CFG.get("max_per_population", 4))
SNAPP_MIN_SAMPLES_LOCUS = int(SNAPP_CFG.get("min_samples_locus", 4))
SNAPP_MISSINGNESS_THRESHOLD = (
    float(SNAPP_CFG.get("missingness_threshold"))
    if SNAPP_CFG.get("missingness_threshold") is not None
    else None
)
_snapp_exclude_cfg = SNAPP_CFG.get("exclude_samples") or []
if isinstance(_snapp_exclude_cfg, str):
    SNAPP_EXCLUDE_SAMPLES = [s.strip() for s in _snapp_exclude_cfg.split(",") if s.strip()]
else:
    SNAPP_EXCLUDE_SAMPLES = [str(s).strip() for s in _snapp_exclude_cfg if str(s).strip()]
SNAPP_POPULATION_ALIASES = {
    str(k): str(v)
    for k, v in (SNAPP_CFG.get("population_aliases", {}) or {}).items()
}
SNAPP_OUTGROUP_ALIASES = {
    str(k): str(v)
    for k, v in (SNAPP_CFG.get("outgroup_aliases", {}) or {}).items()
}
SNAPP_INCLUDE_OUTGROUPS = bool(SNAPP_CFG.get("include_outgroups", True))
SNAPP_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=SNAPP_INCLUDE_OUTGROUPS,
    selected_ids=SNAPP_CFG.get("outgroup_samples"),
)

RAXML_CFG = _config_section("raxml")
RAXML_OUTGROUP_SPECIES = _parse_species_list(RAXML_CFG.get("outgroup_species"))

ANGSD_GLOBAL_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=ANGSD_GLOBAL_INCLUDE_OUTGROUPS,
)
ANGSD_RAXML_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=ANGSD_RAXML_INCLUDE_OUTGROUPS,
    species=RAXML_OUTGROUP_SPECIES,
)
TREEMIX_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=TREEMIX_INCLUDE_OUTGROUPS,
    exclude_ids=TREEMIX_EXCLUDE_OUTGROUPS,
    species=TREEMIX_OUTGROUP_SPECIES,
)
NGSDIST_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=NGSDIST_INCLUDE_OUTGROUPS,
)

ACTIVE_OUTGROUP_SAMPLE_IDS = sorted(
    set(
        ANGSD_GLOBAL_OUTGROUP_IDS
        + ANGSD_RAXML_OUTGROUP_IDS
        + TREEMIX_OUTGROUP_IDS
        + NGSDIST_OUTGROUP_IDS
        + ABBABABA2_OUTGROUP_IDS
        + SNAPP_OUTGROUP_IDS
    )
)

SNAPP_CONSTRAINTS_CFG = (SNAPP_CFG.get("constraints", {}) or {})
SNAPP_CONSTRAINT_TYPE = SNAPP_CONSTRAINTS_CFG.get("type") or SNAPP_CONSTRAINTS_CFG.get("placement")
_taxa_cfg = SNAPP_CONSTRAINTS_CFG.get("taxa")
if isinstance(_taxa_cfg, (list, tuple)):
    SNAPP_CONSTRAINT_TAXA = ",".join(str(x) for x in _taxa_cfg)
else:
    SNAPP_CONSTRAINT_TAXA = _taxa_cfg
SNAPP_RUN_MAP = (SNAPP_CONSTRAINTS_CFG.get("runs", {}) or {})
SNAPP_RUN_IDS = sorted(list(SNAPP_RUN_MAP.keys()))

SNAPP_PREP_CFG = (SNAPP_CFG.get("snapp_prep", {}) or {})
SNAPP_MCMC_LENGTH = int(SNAPP_PREP_CFG.get("mcmc_length", 500000))
SNAPP_TOPOLOGY_WEIGHT = float(SNAPP_PREP_CFG.get("topology_weight", 1.0))
SNAPP_PREP_EXTRA = SNAPP_PREP_CFG.get("extra_args", "").strip()

SNAPP_LOG_PREFIX = SNAPP_CFG.get("log_prefix", "snapp")

#========================#
#### ANALYSIS STATES  ####
#========================#
ANGSD_INTERSECT_ENABLED = _is_enabled(ANGSD_INTERSECT_CFG, True)
ANGSD_GLOBAL_ENABLED = _is_enabled(ANGSD_GLOBAL_CFG, True)
SFS_ENABLED = _is_enabled(SFS_CFG, True)
NGSRELATE_ENABLED = _is_enabled(_config_section("ngsrelate"), True)
NGSLD_ENABLED = _is_enabled(_config_section("ngsld"), True)
PCANGSD_ENABLED = _is_enabled(_config_section("pcangsd"), True)
NGSADMIX_ENABLED = _is_enabled(_config_section("ngsadmix"), True)
NGSDIST_ENABLED = _is_enabled(NGSDIST_CFG, True)
RAXML_ENABLED = _is_enabled(_config_section("raxml"), True)
TREEMIX_ENABLED = _is_enabled(TREEMIX_CFG, True)
SLICE_OUTGROUPS_ENABLED = _is_enabled(_config_section("slice_outgroups"), True)

STRUCTURE_ANALYSES_ENABLED = PCANGSD_ENABLED or NGSADMIX_ENABLED
UNRELATED_ANALYSES_ENABLED = STRUCTURE_ANALYSES_ENABLED or RAXML_ENABLED or TREEMIX_ENABLED
OUTGROUP_ANALYSES_ENABLED = (
    (ANGSD_GLOBAL_ENABLED and bool(ANGSD_GLOBAL_OUTGROUP_IDS))
    or (RAXML_ENABLED and bool(ANGSD_RAXML_OUTGROUP_IDS))
    or (TREEMIX_ENABLED and bool(TREEMIX_OUTGROUP_IDS))
    or (NGSDIST_ENABLED and bool(NGSDIST_OUTGROUP_IDS))
    or (SNAPP_ENABLED and bool(SNAPP_OUTGROUP_IDS))
    or (ABBABABA2_ENABLED and bool(ABBABABA2_OUTGROUP_IDS))
)

ANGSD_INTERSECT_ACTIVE = ANGSD_INTERSECT_ENABLED or ANGSD_GLOBAL_ENABLED or SFS_ENABLED or RAXML_ENABLED or TREEMIX_ENABLED or ABBABABA2_ENABLED or SNAPP_ENABLED or NGSDIST_ENABLED
ANGSD_GLOBAL_ACTIVE = ANGSD_GLOBAL_ENABLED or NGSRELATE_ENABLED or NGSLD_ENABLED or STRUCTURE_ANALYSES_ENABLED or SNAPP_ENABLED
NGSRELATE_ACTIVE = NGSRELATE_ENABLED or UNRELATED_ANALYSES_ENABLED
NGSLD_ACTIVE = NGSLD_ENABLED or TREEMIX_ENABLED or (SFS_ENABLED and SFS_NEEDS_UNLINKED_SITES)
GLOBAL_UNRELATED_UNLINKED_ACTIVE = STRUCTURE_ANALYSES_ENABLED
SLICE_OUTGROUPS_ACTIVE = SLICE_OUTGROUPS_ENABLED and bool(ACTIVE_OUTGROUP_SAMPLE_IDS)
