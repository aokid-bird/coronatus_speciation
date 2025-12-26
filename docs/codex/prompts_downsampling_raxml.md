Because raxml-ng takes really long if many samples are involved, for a test run and other purposes, the samples needed to be downsampled for same cases. In the workflow/rules/angsd_raxml.smk, can you implement this? Please consider the following.

1) Users can specify the maximum number of samples per defined populations (in the config) to be contained in the config file.
2) Users can specify samples that they do not want to include in the downsampled subsets.
3) Randomly choose the samples up to the maximum number defined (1). If sample size is fewer than the max, use all.
4) Let users to use all the samples by specifying in the config file. If this variable is specified, then the max number will be ignored.
5) Please stick with the present raxml-ng run flows. For the downsampling and other settings I just listed here, please implement them with minimum changes, especially in the catg file reconstruction and raxml-ng runs. Maybe (but please think well on your own), implement them in the upstream of the rules, before angsd_raxml.
6) Please update this document to include summaries about what you've done as the response to my prompts and notes anything you should tell users about the usage of the functionality.   

---

## Codex response (2024-05-26)

- Added `rule make_bamlist_raxml_downsampled` so RAxML analyses work from a per-population downsampled bamlist while keeping the existing `angsd_raxml`/CATG/RAxML-ng flow unchanged.
- Extended `workflow/rules/common.smk` to expose downsampling parameters from `angsd_raxml.downsampling` and pass them into the new rule.
- `workflow/rules/angsd_raxml.smk` now builds the analysis bamlist from the downsampled ingroup set before appending optional outgroups.

### Usage notes

- Configure downsampling in the active YAML config under `angsd_raxml.downsampling`. Supported keys:
  - `max_per_population`: integer applied to every population or a mapping like `{hardwickii: 6, megala_e: 4, default: 5}`. Omit or set `null` to use all available samples.
  - `exclude_samples`: list (or comma-delimited string) of sample IDs to always drop before sampling.
  - `use_all_samples`: set `true` to bypass the per-population max while still respecting `exclude_samples`.
  - `seed`: optional integer to control the deterministic random draw (defaults to `0` when omitted).
- The resulting bamlist lives at `results/bamlists/<output_prefix>/raxml_downsampled/bamlist.txt`; `raxml_analysis/bamlist.txt` continues to be the input for ANGSD/RAxML, now including any requested outgroups plus the downsampled ingroup.
