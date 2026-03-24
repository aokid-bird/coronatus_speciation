"""
Shared config normalization and lightweight global settings.

This file is responsible for:
- backfilling legacy flat config keys from the grouped schema
- exposing a small set of global workflow settings
- parsing generic analysis toggles and SFS settings
"""

import os
import re
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import pandas as pd

_CFG_MISSING = object()


# Generic nested-config helpers used throughout the shared layer.
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


# Backfill old top-level keys so existing rule files can keep working while
# config defaults are organized under grouped sections.
for target_keys, source_keys, default in [
    (("environment",), ("project", "environment"), "cluster"),
    (("threads",), ("project", "threads"), 1),
    (("output_prefix",), ("project", "output_prefix"), "defaults"),
    (("populations",), ("project", "populations"), []),
    (("group_col",), ("project", "group_col"), None),
    (("population_labels",), ("project", "population_labels"), None),
    (("bam_dir",), ("paths", "bam_dir"), "results/bwa"),
    (("singularity_dir",), ("paths", "singularity_dir"), "${HOME}/envs/singularity"),
    (("transfer",), ("paths", "transfer"), {}),
    (("samples",), ("inputs", "samples_tsv"), "data/samples.tsv"),
    (("outgroups",), ("inputs", "outgroups_tsv"), "data/outgroup.tsv"),
    (("references_tsv",), ("inputs", "references_tsv"), "data/references.tsv"),
    (("reference_download_method",), ("reference", "download_method"), "datasets"),
    (("reference_dir",), ("reference", "storage_dir"), "data/reference"),
    (("outgroups_seq_column",), ("metadata", "outgroups_seq_column"), "sequencer"),
    (("longread_keywords",), ("metadata", "longread_keywords"), ["PacBio", "ONT", "Nanopore"]),
    (("angsd_common_args",), ("angsd", "common_args"), ""),
    (("angsd_args",), ("angsd", "args"), {}),
    (("minIndRatio",), ("angsd", "min_ind_ratio"), {}),
]:
    _backfill_config(target_keys, source_keys, default=default)

# Backfill grouped analysis configs to their legacy top-level names.
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

output_prefix = config["output_prefix"]
groups = config["populations"]
group_col = config["group_col"]
LEGACY_GLOBAL_THREADS = int(config.get("threads", 1))
ENVIRONMENT = config.get("environment", "cluster")


# Small shared accessors used by the rest of the workflow.
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


def _is_enabled(section, default=True):
    return bool(section.get("enabled", default))


def _parse_list(raw):
    if raw is None:
        return []
    if isinstance(raw, str):
        return [item.strip() for item in raw.split(",") if item.strip()]
    return [str(item).strip() for item in raw if str(item).strip()]


def _parse_optional_list(raw):
    if raw is None:
        return None
    parsed = _parse_list(raw)
    return parsed if parsed else []


def _parse_species_list(raw):
    return _parse_list(raw)


def _normalise_seed(value):
    if value is None or (isinstance(value, str) and value.strip().lower() in {"", "none"}):
        return None
    return int(value)


_POP_LABEL_CFG = config.get("population_labels") or {}
if isinstance(_POP_LABEL_CFG, dict):
    POPULATION_LABELS = {str(key): str(value) for key, value in _POP_LABEL_CFG.items()}
else:
    POPULATION_LABELS = {}


def get_population_label(pop):
    return POPULATION_LABELS.get(str(pop), str(pop))


SFS_CFG = config.get("sfs_analysis", {}) or {}


# Normalize simple token-like config values used in SFS target naming.
def _norm_token(value, default):
    token = str(value).strip().lower().replace(" ", "_") if value is not None else None
    return token if token else default


SFS_SITE_FILTERS = [
    _norm_token(value, "linked")
    for value in SFS_CFG.get("site_filters", ["linked"]) or ["linked"]
]
SFS_FOLD_STATES = [
    _norm_token(value, "fold")
    for value in SFS_CFG.get("fold_states", ["fold"]) or ["fold"]
]
SFS_INCLUDE_MONOMORPHIC = bool(SFS_CFG.get("include_monomorphic", True))
_SFS_REAL_CFG = SFS_CFG.get("realSFS", {}) or {}
SFS_MAXITER = int(_SFS_REAL_CFG.get("maxiter", 50000))
SFS_TOLE = str(_SFS_REAL_CFG.get("tole", "1e-6"))
SFS_NEEDS_UNLINKED_SITES = "unlinked" in {value.lower() for value in SFS_SITE_FILTERS}
SFS_PAIRWISE_COMBOS = list(combinations(groups, 2)) if len(groups) >= 2 else []
