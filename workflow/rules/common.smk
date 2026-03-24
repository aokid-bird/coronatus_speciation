"""
Compatibility shim for the refactored shared workflow includes.

Prefer including the narrower common_*.smk files directly from the main Snakefile.
"""

include: "rules/common_config.smk"
include: "rules/common_paths.smk"
include: "rules/common_samples.smk"
include: "rules/common_analysis.smk"
include: "rules/common_rules.smk"
