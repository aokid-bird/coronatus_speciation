import pandas as pd
from pathlib import Path
from itertools import combinations
import re
import os

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
config_bam_base = config['bam_dir']
config_bam_dir = f"{config_bam_base}/{output_prefix}"
OUTGROUP_SLICED_DIR = f"{config_bam_dir}/outgroups_sliced"
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
REFERENCE_DIR = config.get("reference_dir", "data/reference")
REFERENCE_METHOD = config.get("reference_download_method", "datasets")

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
INGROUP_READS_DIR = (config.get("reads", {}) or {}).get("ingroup_dir", "data/raw/ingroup")

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

def _parse_exclude_list(raw):
    if raw is None:
        return []
    if isinstance(raw, str):
        return [s.strip() for s in raw.split(",") if s.strip()]
    return [str(s).strip() for s in raw if str(s).strip()]

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

_abb_sel = ABBABABA2_CFG.get("outgroup_samples")
if _abb_sel is None:
    ABBABABA2_OUTGROUP_IDS = list(OUTGROUP_SAMPLE_IDS)
elif isinstance(_abb_sel, str):
    ABBABABA2_OUTGROUP_IDS = [s.strip() for s in _abb_sel.split(",") if s.strip()]
else:
    ABBABABA2_OUTGROUP_IDS = [str(s).strip() for s in _abb_sel if str(s).strip()]
ABBABABA2_OUTGROUP_IDS = [sid for sid in ABBABABA2_OUTGROUP_IDS if sid in OUTGROUP_SAMPLE_IDS]

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

ABBABABA2_KHROMA_LIB = str(ABBABABA2_CFG.get("khroma_lib", "") or "")

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
_snapp_outgroup_sel = SNAPP_CFG.get("outgroup_samples")
if SNAPP_INCLUDE_OUTGROUPS:
    if isinstance(_snapp_outgroup_sel, str):
        SNAPP_OUTGROUP_IDS = [s.strip() for s in _snapp_outgroup_sel.split(",") if s.strip()]
    elif isinstance(_snapp_outgroup_sel, (list, tuple)):
        SNAPP_OUTGROUP_IDS = [str(s) for s in _snapp_outgroup_sel]
    else:
        SNAPP_OUTGROUP_IDS = list(OUTGROUP_SAMPLE_IDS)
    SNAPP_OUTGROUP_IDS = [sid for sid in SNAPP_OUTGROUP_IDS if sid in OUTGROUP_SAMPLE_IDS]
else:
    SNAPP_OUTGROUP_IDS = []

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
