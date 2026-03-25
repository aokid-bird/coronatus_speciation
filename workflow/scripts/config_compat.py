from __future__ import annotations

from pathlib import Path

import yaml

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


def _backfill_config(config, target_keys, *source_options, default=_CFG_MISSING):
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


def apply_external_config_compat(config: dict) -> dict:
    """Backfill legacy flat keys for standalone helpers outside Snakemake."""
    for target_keys, source_keys, default in [
        (("environment",), ("project", "environment"), "cluster"),
        (("threads",), ("project", "threads"), 1),
        (("output_prefix",), ("project", "output_prefix"), "defaults"),
        (("populations",), ("project", "populations"), []),
        (("group_col",), ("project", "group_col"), None),
        (("population_labels",), ("project", "population_labels"), None),
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
        _backfill_config(config, target_keys, source_keys, default=default)

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
        _backfill_config(config, (analysis_name,), ("analyses", analysis_name), default={})

    _backfill_config(config, ("sfs_analysis",), ("analyses", "sfs"), default={})
    return config


def load_config_with_compat(config_path: str | Path) -> dict:
    with open(config_path) as handle:
        config = yaml.safe_load(handle) or {}
    return apply_external_config_compat(config)
