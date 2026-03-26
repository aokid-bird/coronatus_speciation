"""
Shared filesystem, storage, reference, and mapper path setup.

This file keeps path construction in one place so rule files can refer to
resolved directories and reference assets without repeating config logic.
"""

STORAGE_CFG = config.get("storage", {}) or {}
TRANSFER_CFG = config.get("transfer", {}) or {}


# Read nested storage paths with sensible defaults.
def _storage_path(*keys, default=None):
    value = STORAGE_CFG
    for key in keys:
        if not isinstance(value, dict):
            return default
        value = value.get(key)
    return default if value in (None, "") else value


def _storage_root_for(environment):
    value = _config_value(STORAGE_CFG, "roots", environment, default=_CFG_MISSING)
    if value is _CFG_MISSING or value in (None, ""):
        return None
    return str(Path(os.path.expandvars(str(value))).expanduser())


def _storage_shared_path(*keys):
    value = _config_value(STORAGE_CFG, "shared", *keys, default=_CFG_MISSING)
    if value is _CFG_MISSING or value in (None, ""):
        return None
    return str(value)


def _resolve_storage_path(explicit_keys, shared_keys=None, default=None):
    explicit = _storage_path(*explicit_keys, default=None)
    if explicit not in (None, ""):
        return explicit
    if ACTIVE_STORAGE_ROOT and shared_keys:
        shared = _storage_shared_path(*shared_keys)
        if shared:
            return str(Path(ACTIVE_STORAGE_ROOT) / shared)
    return default


def _translate_between_storage_roots(path, source_root, target_root):
    if not source_root or not target_root:
        return None
    candidate = Path(str(path)).expanduser()
    try:
        relative = candidate.relative_to(Path(source_root))
    except ValueError:
        return None
    return str(Path(target_root) / relative)


def rebase_path_to_active_storage_root(path):
    if not path:
        return path
    candidate = Path(str(path)).expanduser()
    if not candidate.is_absolute():
        return str(path)
    translated = _translate_between_storage_roots(candidate, LOCAL_STORAGE_ROOT, ACTIVE_STORAGE_ROOT)
    return translated or str(candidate)


INGROUP_FASTQ_BASENAME_COL = str(
    _config_value(config, "reads", "ingroup_metadata", "basename_col", default="sample")
)
INGROUP_BAM_ALIAS_ACTIVE = INGROUP_FASTQ_BASENAME_COL != "sample"
LOCAL_STORAGE_ROOT = _storage_root_for("local")
CLUSTER_STORAGE_ROOT = _storage_root_for("cluster")
ACTIVE_STORAGE_ROOT = _storage_root_for(ENVIRONMENT)


def _legacy_bam_base():
    value = _config_value(config, "paths", "bam_dir", default=_CFG_MISSING)
    if value is _CFG_MISSING or value in (None, ""):
        return "results/bwa"
    return str(value)


LEGACY_BAM_BASE = _legacy_bam_base()
LEGACY_INGROUP_BAM_DIR = (
    f"{LEGACY_BAM_BASE}/{output_prefix}/preprocess"
    if INGROUP_BAM_ALIAS_ACTIVE
    else f"{LEGACY_BAM_BASE}/{output_prefix}"
)
INGROUP_BAM_STORAGE_DIR = _resolve_storage_path(
    ("bam", "ingroup_dir"),
    ("bam", "ingroup_dir"),
    default=LEGACY_INGROUP_BAM_DIR,
)
config_bam_dir = (
    f"{LEGACY_BAM_BASE}/{output_prefix}"
    if INGROUP_BAM_ALIAS_ACTIVE
    else INGROUP_BAM_STORAGE_DIR
)
OUTGROUP_BAM_DIR = _resolve_storage_path(
    ("bam", "outgroup_dir"),
    ("bam", "outgroup_dir"),
    default=(
        f"{LEGACY_BAM_BASE}/{output_prefix}/outgroups"
        if _config_value(config, "paths", "bam_dir", default=_CFG_MISSING) is not _CFG_MISSING
        else f"{INGROUP_BAM_STORAGE_DIR}/outgroups"
    ),
)
OUTGROUP_SLICED_DIR = f"results/outgroups_sliced/{output_prefix}"

CLUSTER_STORAGE_ROOT_DEFAULT = str(
    TRANSFER_CFG.get("cluster_storage_root")
    or CLUSTER_STORAGE_ROOT
    or "/lfs/aokid"
)
PROJECT_ROOT = Path.cwd().resolve()
config_singularity_dir = os.path.expandvars(config["singularity_dir"])

REF_CONFIG_PATH = config.get("reference", {}).get("fasta", config.get("ref"))
REFERENCES_TSV = config.get("references_tsv", "data/references.tsv")
REFERENCE_DIR = _resolve_storage_path(
    ("reference", "dir"),
    ("reference", "dir"),
    default=config.get("reference_dir", "data/reference"),
)
REFERENCE_METHOD = config.get("reference_download_method", "datasets")
OUTGROUP_RAW_DIR = _resolve_storage_path(
    ("outgroup", "raw_dir"),
    ("outgroup", "raw_dir"),
    default="data/raw/outgroup",
)
OUTGROUP_MERGED_DIR = _resolve_storage_path(
    ("outgroup", "merged_dir"),
    ("outgroup", "merged_dir"),
    default="data/merged/outgroup",
)
INGROUP_TRIM_DIR = _resolve_storage_path(
    ("derived", "ingroup", "trim_dir"),
    ("derived", "ingroup", "trim_dir"),
    default="results/trimmomatic/ingroup",
)
OUTGROUP_TRIM_DIR = _resolve_storage_path(
    ("derived", "outgroup", "trim_dir"),
    ("derived", "outgroup", "trim_dir"),
    default="results/trimmomatic/outgroup",
)
INGROUP_QC_BASE_DIR = _resolve_storage_path(
    ("derived", "ingroup", "qc_dir"),
    ("derived", "ingroup", "qc_dir"),
    default="results/qc/ingroup",
)
OUTGROUP_QC_BASE_DIR = _resolve_storage_path(
    ("derived", "outgroup", "qc_dir"),
    ("derived", "outgroup", "qc_dir"),
    default="results/qc/outgroup",
)
INGROUP_MAP_TMP_DIR = _resolve_storage_path(
    ("derived", "ingroup", "mapping_tmp_dir"),
    ("derived", "ingroup", "mapping_tmp_dir"),
    default=f"results/mapping/{output_prefix}/ingroup",
)
OUTGROUP_MAP_TMP_DIR = _resolve_storage_path(
    ("derived", "outgroup", "mapping_tmp_dir"),
    ("derived", "outgroup", "mapping_tmp_dir"),
    default=f"results/mapping/{output_prefix}/outgroup",
)
OUTGROUP_LR_FILTER_DIR = _resolve_storage_path(
    ("derived", "outgroup", "longread_filter_dir"),
    ("derived", "outgroup", "longread_filter_dir"),
    default="results/longread/filter",
)


def resolve_ingroup_reads_dir():
    explicit = _config_value(config, "reads", "ingroup_dir", default=_CFG_MISSING)
    if explicit is not _CFG_MISSING and explicit not in (None, ""):
        return explicit
    shared = _storage_shared_path("reads", "ingroup_dir")
    if ACTIVE_STORAGE_ROOT and shared:
        return str(Path(ACTIVE_STORAGE_ROOT) / shared)
    return "data/raw/ingroup"


# Standard storage-manifest filenames for managed directories.
def manifest_paths(directory):
    return (
        str(Path(directory) / "README.md"),
        str(Path(directory) / "provenance.yaml"),
    )


def _is_directory_target_writable(directory):
    probe = Path(directory)
    while not probe.exists():
        if probe.parent == probe:
            return False
        probe = probe.parent
    return probe.is_dir() and os.access(probe, os.W_OK)


def manifest_directory(directory, fallback_directory=None):
    if _is_directory_target_writable(directory):
        return directory
    if fallback_directory:
        return fallback_directory
    return directory


REFERENCE_MANIFEST_DIR = manifest_directory(
    REFERENCE_DIR,
    fallback_directory=f"results/reference/{output_prefix}/source_metadata",
)


# Write a small README and provenance YAML for workflow-managed storage paths.
def write_storage_manifest(directory, title, producer, details=None):
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


# Resolve reference sidecar index names for compressed and uncompressed FASTA.
def derive_indexes(fasta):
    fasta_path = Path(fasta)
    if fasta_path.suffix == ".gz":
        return f"{fasta_path}.fai", f"{fasta_path}.gzi"
    return f"{fasta_path}.fai", None


REF_FROM_METADATA = False
REFERENCE_ACCESSION = None
REFERENCE_NAME = None
REFERENCE_URL = ""
REFERENCE_ZIP = None

if REF_CONFIG_PATH:
    REF = REF_CONFIG_PATH
else:
    if ENVIRONMENT != "local":
        raise ValueError(
            "On cluster environment, set reference.fasta in config to the transferred reference."
        )
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

FAI_DEFAULT, GZI = derive_indexes(REF)
FAI = config.get("reference", {}).get("fai") or FAI_DEFAULT
REF_CHR = (
    f"{REFERENCE_DIR}/{Path(REF).stem}.chr"
    if REF_FROM_METADATA
    else f"{Path(REF).with_suffix('').with_suffix('').as_posix()}.chr"
)


# Path utilities used for transfer-script generation.
def _path_is_within(path, base):
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False


def is_external_storage_path(path):
    candidate = Path(str(path)).expanduser()
    if not candidate.is_absolute():
        return False
    return not _path_is_within(candidate, PROJECT_ROOT)


def remote_mirror_path(path, cluster_root=None):
    target_root = cluster_root or CLUSTER_STORAGE_ROOT_DEFAULT
    return str(Path(target_root) / remote_mirror_subpath(path))


def remote_mirror_subpath(path):
    candidate = Path(str(path)).expanduser()
    for source_root in (LOCAL_STORAGE_ROOT, ACTIVE_STORAGE_ROOT):
        if not source_root:
            continue
        try:
            return str(candidate.relative_to(Path(source_root)))
        except ValueError:
            continue
    return str(candidate).lstrip("/")


def local_mirror_path(path, local_root=None):
    target_root = local_root or LOCAL_STORAGE_ROOT
    subpath = remote_mirror_subpath(path)
    if target_root:
        return str(Path(target_root) / subpath)
    return str(Path(str(path)).expanduser())


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
        INGROUP_BAM_STORAGE_DIR,
        OUTGROUP_BAM_DIR,
    ]
    return sorted({path for path in candidates if is_external_storage_path(path)})


# Build per-scope mapper config by overlaying ingroup/outgroup settings on the
# shared mapping section.
MAPPING_CFG = config.get("mapping", {}) or {}
MAPPING_INGROUP_CFG = {k: v for k, v in MAPPING_CFG.items() if k not in {"ingroup", "outgroup"}}
MAPPING_OUTGROUP_CFG = dict(MAPPING_INGROUP_CFG)
MAPPING_INGROUP_CFG.update(MAPPING_CFG.get("ingroup", {}) or {})
MAPPING_OUTGROUP_CFG.update(MAPPING_CFG.get("outgroup", {}) or {})
MAPPING_INGROUP_THREADS = _resolve_named_threads(6, "mapping", "ingroup", legacy_section=MAPPING_INGROUP_CFG)
MAPPING_OUTGROUP_THREADS = _resolve_named_threads(6, "mapping", "outgroup", legacy_section=MAPPING_OUTGROUP_CFG)
MAPPER_INGROUP = MAPPING_INGROUP_CFG.get("mapper", config.get("mapper", "bwa"))
MAPPER_OUTGROUP = MAPPING_OUTGROUP_CFG.get("mapper", MAPPER_INGROUP)

ref_path = Path(REF)
REF_UNZIPPED = str(ref_path.with_suffix("")) if ref_path.suffix == ".gz" else str(ref_path)
REF_INDEX_DIR = (config.get("reference", {}) or {}).get("index_dir")
REF_INDEX_PREFIX_BASE = (
    str(Path(REF_INDEX_DIR) / Path(REF_UNZIPPED).name)
    if REF_INDEX_DIR
    else REF_UNZIPPED
)

REF_INDEX_PREFIX_BWA = REF_INDEX_PREFIX_BASE
REF_INDEX_PREFIX_BWAMEM2 = REF_INDEX_PREFIX_BASE
if MAPPER_INGROUP != MAPPER_OUTGROUP:
    REF_INDEX_PREFIX_BWA = f"{REF_INDEX_PREFIX_BASE}.bwa"
    REF_INDEX_PREFIX_BWAMEM2 = f"{REF_INDEX_PREFIX_BASE}.bwamem2"

BWA_INDEX_MAIN = f"{REF_INDEX_PREFIX_BWA}.bwt"
BWAMEM2_INDEX_MAIN = f"{REF_INDEX_PREFIX_BWAMEM2}.0123"


# Select the representative mapper index output used as a dependency target.
def _ref_index_for(mapper):
    return BWAMEM2_INDEX_MAIN if mapper == "bwa-mem2" else BWA_INDEX_MAIN


REF_MAP_INDEX_INGROUP = _ref_index_for(MAPPER_INGROUP)
REF_MAP_INDEX_OUTGROUP = _ref_index_for(MAPPER_OUTGROUP)
REF_MAP_INDEX_MAIN = REF_MAP_INDEX_INGROUP
REF_MAP_INDEXES = sorted({REF_MAP_INDEX_INGROUP, REF_MAP_INDEX_OUTGROUP})
REF_MAP_ARG_BWA = REF_INDEX_PREFIX_BWA
REF_MAP_ARG_BWAMEM2 = REF_INDEX_PREFIX_BWAMEM2
REF_MAP_ARG_INGROUP = REF_MAP_ARG_BWAMEM2 if MAPPER_INGROUP == "bwa-mem2" else REF_MAP_ARG_BWA
REF_MAP_ARG_OUTGROUP = REF_MAP_ARG_BWAMEM2 if MAPPER_OUTGROUP == "bwa-mem2" else REF_MAP_ARG_BWA
REF_MAP_ARG = REF_MAP_ARG_INGROUP
REF_INDEX_PREFIX = REF_INDEX_PREFIX_BASE
MAPPER = MAPPER_INGROUP
