"""
Shared analysis-level config parsing and workflow state.

This file keeps the analysis toggles, thread settings, shared selection state,
and cross-analysis activation logic in one place.
"""

# Load grouped configs for each analysis family once.
ANGSD_INTERSECT_CFG = config.get("angsd_intersect", {}) or {}
ANGSD_GLOBAL_CFG = config.get("angsd_global", {}) or {}
ANGSD_RAXML_CFG = config.get("angsd_raxml", {}) or {}
ANGSD_SFS_CFG = config.get("angsd_sfs", {}) or {}
NGSRELATE_CFG = _config_section("ngsrelate")
NGSLD_CFG = _config_section("ngsld")
PCANGSD_CFG = _config_section("pcangsd")
NGSADMIX_CFG = _config_section("ngsadmix")
NGSDIST_CFG = _config_section("ngsdist")
RAXML_CFG = _config_section("raxml")
TREEMIX_CFG = config.get("treemix", {}) or {}
SLICE_OUTGROUPS_CFG = _config_section("slice_outgroups")
ABBABABA2_CFG = config.get("abbababa2", {}) or {}
SNAPP_CFG = config.get("snapp", {}) or {}

# Shared per-analysis toggles and thread settings referenced by many rule files.
ANGSD_GLOBAL_INCLUDE_OUTGROUPS = bool(ANGSD_GLOBAL_CFG.get("include_outgroups", False))
ANGSD_RAXML_INCLUDE_OUTGROUPS = bool(ANGSD_RAXML_CFG.get("include_outgroups", False))
NGSDIST_INCLUDE_OUTGROUPS = bool(NGSDIST_CFG.get("include_outgroups", ANGSD_GLOBAL_INCLUDE_OUTGROUPS))
TREEMIX_INCLUDE_OUTGROUPS = bool(TREEMIX_CFG.get("include_outgroups", True))
SNAPP_INCLUDE_OUTGROUPS = bool(SNAPP_CFG.get("include_outgroups", True))

ANGSD_INTERSECT_THREADS = _resolve_threads(ANGSD_INTERSECT_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_GLOBAL_THREADS = _resolve_threads(ANGSD_GLOBAL_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_GLOBAL_UNRELATED_THREADS = _resolve_threads(
    _config_section("angsd_global_unrelated_unlinked"),
    ANGSD_GLOBAL_THREADS,
)
ANGSD_SFS_THREADS = _resolve_threads(ANGSD_SFS_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_RAXML_THREADS = _resolve_threads(ANGSD_RAXML_CFG, LEGACY_GLOBAL_THREADS)
ANGSD_SNAPP_THREADS = _resolve_threads(SNAPP_CFG, ANGSD_GLOBAL_THREADS)
SLICE_OUTGROUPS_THREADS = _resolve_threads(SLICE_OUTGROUPS_CFG, 6)
NGSLD_THREADS = _resolve_threads(NGSLD_CFG, LEGACY_GLOBAL_THREADS)
REALSFS_THREADS = _resolve_threads(_config_section("sfs_analysis", "realSFS"), 10, legacy_fallback=False)
QC_INGROUP_THREADS = _resolve_threads(_config_section("qc", "ingroup", "trimmomatic"), 4)
QC_OUTGROUP_THREADS = _resolve_threads(_config_section("qc", "outgroup", "trimmomatic"), 4)

kin_thr = NGSRELATE_CFG["kinship_threshold"]

TREEMIX_MODE = TREEMIX_CFG.get("mode", "pop")
TREEMIX_ROOT_LABEL = TREEMIX_CFG.get("root_label")
TREEMIX_EXCLUDE_SAMPLES = _parse_list(TREEMIX_CFG.get("exclude_samples"))
TREEMIX_OUTGROUP_SPECIES = _parse_species_list(TREEMIX_CFG.get("outgroup_species"))

ANGSD_RAXML_DOWNSAMPLE_CFG = ANGSD_RAXML_CFG.get("downsampling", {}) or {}


# Normalize downsampling limits, allowing ints, dicts, and null-like values.
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


ANGSD_RAXML_DOWNSAMPLE_MAX = None
try:
    ANGSD_RAXML_DOWNSAMPLE_MAX = _parse_max_per_population(
        ANGSD_RAXML_DOWNSAMPLE_CFG.get("max_per_population")
    )
except (TypeError, ValueError):
    ANGSD_RAXML_DOWNSAMPLE_MAX = None

ANGSD_RAXML_DOWNSAMPLE_EXCLUDE = _parse_list(
    ANGSD_RAXML_DOWNSAMPLE_CFG.get("exclude_samples")
)
ANGSD_RAXML_DOWNSAMPLE_USE_ALL = bool(
    ANGSD_RAXML_DOWNSAMPLE_CFG.get("use_all_samples", False)
)
try:
    ANGSD_RAXML_DOWNSAMPLE_SEED = _normalise_seed(ANGSD_RAXML_DOWNSAMPLE_CFG.get("seed"))
except (TypeError, ValueError):
    ANGSD_RAXML_DOWNSAMPLE_SEED = None


# Resolve ANGSD minInd ratios with the same fallback order across analyses.
def get_minInd_ratio(key, default=None):
    ratios = config.get("minIndRatio", {}) or {}
    if key in ratios:
        return float(ratios[key])
    if key != "global" and "global" in ratios:
        return float(ratios["global"])
    if key != "intersect" and "intersect" in ratios:
        return float(ratios["intersect"])
    return default


# Shared outgroup selections used by multiple downstream analyses.
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
    exclude_ids=TREEMIX_CFG.get("exclude_outgroups"),
    species=TREEMIX_OUTGROUP_SPECIES,
)
NGSDIST_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=NGSDIST_INCLUDE_OUTGROUPS,
)

ABBABABA2_ENABLED = bool(ABBABABA2_CFG.get("enabled", False))
ABBABABA2_THREADS = int(ABBABABA2_CFG.get("threads", LEGACY_GLOBAL_THREADS))
ABBABABA2_ANGSD_ARGS = str(ABBABABA2_CFG.get("angsd_args", "")).strip()
if not ABBABABA2_ANGSD_ARGS:
    ABBABABA2_ANGSD_ARGS = "-doAbbababa2 1 -doCounts 1 -minMapQ 30 -minQ 20 -baq 2 -useLast 1"
ABBABABA2_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=True,
    selected_ids=ABBABABA2_CFG.get("outgroup_samples"),
)
ABBABABA2_EXCLUDE_SAMPLES = _parse_list(ABBABABA2_CFG.get("exclude_samples"))
_ABBABABA2_LABEL_CFG = ABBABABA2_CFG.get("outgroup_label")
ABBABABA2_OUTGROUP_LABEL_DEFAULT = None
if isinstance(_ABBABABA2_LABEL_CFG, dict):
    ABBABABA2_OUTGROUP_LABELS = {
        str(key): str(value)
        for key, value in _ABBABABA2_LABEL_CFG.items()
        if str(value).strip()
    }
else:
    ABBABABA2_OUTGROUP_LABELS = {}
    if _ABBABABA2_LABEL_CFG is not None:
        label = str(_ABBABABA2_LABEL_CFG).strip()
        if label:
            ABBABABA2_OUTGROUP_LABEL_DEFAULT = label


# Resolve the label that ABBABABA2 should use for an outgroup sample.
def resolve_abbababa2_outgroup_label(sample_id):
    sid = str(sample_id)
    if sid in ABBABABA2_OUTGROUP_LABELS:
        return ABBABABA2_OUTGROUP_LABELS[sid]
    if ABBABABA2_OUTGROUP_LABEL_DEFAULT:
        return ABBABABA2_OUTGROUP_LABEL_DEFAULT
    return OUTGROUP_TAXON_LABELS.get(sid) or sid


SNAPP_ENABLED = bool(SNAPP_CFG.get("enabled", False))
SNAPP_MAX_PER_POP = int(SNAPP_CFG.get("max_per_population", 4))
SNAPP_MIN_SAMPLES_LOCUS = int(SNAPP_CFG.get("min_samples_locus", 4))
SNAPP_MISSINGNESS_THRESHOLD = (
    float(SNAPP_CFG.get("missingness_threshold"))
    if SNAPP_CFG.get("missingness_threshold") is not None
    else None
)
SNAPP_EXCLUDE_SAMPLES = _parse_list(SNAPP_CFG.get("exclude_samples"))
SNAPP_POPULATION_ALIASES = {
    str(key): str(value)
    for key, value in (SNAPP_CFG.get("population_aliases", {}) or {}).items()
}
SNAPP_OUTGROUP_ALIASES = {
    str(key): str(value)
    for key, value in (SNAPP_CFG.get("outgroup_aliases", {}) or {}).items()
}
SNAPP_OUTGROUP_IDS = analysis_outgroup_sample_ids(
    include_outgroups=SNAPP_INCLUDE_OUTGROUPS,
    selected_ids=SNAPP_CFG.get("outgroup_samples"),
)
SNAPP_CONSTRAINTS_CFG = SNAPP_CFG.get("constraints", {}) or {}
SNAPP_CONSTRAINT_TYPE = SNAPP_CONSTRAINTS_CFG.get("type") or SNAPP_CONSTRAINTS_CFG.get("placement")
_SNAPP_TAXA_CFG = SNAPP_CONSTRAINTS_CFG.get("taxa")
SNAPP_CONSTRAINT_TAXA = (
    ",".join(str(value) for value in _SNAPP_TAXA_CFG)
    if isinstance(_SNAPP_TAXA_CFG, (list, tuple))
    else _SNAPP_TAXA_CFG
)
SNAPP_RUN_MAP = SNAPP_CONSTRAINTS_CFG.get("runs", {}) or {}
SNAPP_RUN_IDS = sorted(SNAPP_RUN_MAP.keys())
SNAPP_PREP_CFG = SNAPP_CFG.get("snapp_prep", {}) or {}
SNAPP_MCMC_LENGTH = int(SNAPP_PREP_CFG.get("mcmc_length", 500000))
SNAPP_TOPOLOGY_WEIGHT = float(SNAPP_PREP_CFG.get("topology_weight", 1.0))
SNAPP_PREP_EXTRA = SNAPP_PREP_CFG.get("extra_args", "").strip()
SNAPP_LOG_PREFIX = SNAPP_CFG.get("log_prefix", "snapp")

# Union of sliced outgroups actually needed by enabled analyses.
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

ANGSD_INTERSECT_ENABLED = _is_enabled(ANGSD_INTERSECT_CFG, True)
ANGSD_GLOBAL_ENABLED = _is_enabled(ANGSD_GLOBAL_CFG, True)
SFS_ENABLED = _is_enabled(SFS_CFG, True)
NGSRELATE_ENABLED = _is_enabled(NGSRELATE_CFG, True)
NGSLD_ENABLED = _is_enabled(NGSLD_CFG, True)
PCANGSD_ENABLED = _is_enabled(PCANGSD_CFG, True)
NGSADMIX_ENABLED = _is_enabled(NGSADMIX_CFG, True)
NGSDIST_ENABLED = _is_enabled(NGSDIST_CFG, True)
RAXML_ENABLED = _is_enabled(RAXML_CFG, True)
TREEMIX_ENABLED = _is_enabled(TREEMIX_CFG, True)
SLICE_OUTGROUPS_ENABLED = _is_enabled(SLICE_OUTGROUPS_CFG, True)

# Derived workflow states used by Snakefile and dependency wiring.
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

ANGSD_INTERSECT_ACTIVE = (
    ANGSD_INTERSECT_ENABLED
    or ANGSD_GLOBAL_ENABLED
    or SFS_ENABLED
    or RAXML_ENABLED
    or TREEMIX_ENABLED
    or ABBABABA2_ENABLED
    or SNAPP_ENABLED
    or NGSDIST_ENABLED
)
ANGSD_GLOBAL_ACTIVE = (
    ANGSD_GLOBAL_ENABLED
    or NGSRELATE_ENABLED
    or NGSLD_ENABLED
    or STRUCTURE_ANALYSES_ENABLED
    or SNAPP_ENABLED
)
NGSRELATE_ACTIVE = NGSRELATE_ENABLED or UNRELATED_ANALYSES_ENABLED
NGSLD_ACTIVE = NGSLD_ENABLED or TREEMIX_ENABLED or (SFS_ENABLED and SFS_NEEDS_UNLINKED_SITES)
GLOBAL_UNRELATED_UNLINKED_ACTIVE = STRUCTURE_ANALYSES_ENABLED
SLICE_OUTGROUPS_ACTIVE = SLICE_OUTGROUPS_ENABLED and bool(ACTIVE_OUTGROUP_SAMPLE_IDS)
