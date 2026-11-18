# CARTaGENE ancestry-stratified GRCh38 LD

This pipeline turns chromosome-level imputed PGEN files into ancestry-specific,
QC-filtered BED files and signed LD correlation matrices using magenpy. Genotypes
stay on the SLURM cluster. Autosomes 1–22 are supported. The implementation follows
[HGDP+1KG](https://github.com/shz9/1KG_HGDP_ld), with early panel extraction to
handle CARTaGENE's approximately 455 million input variants and 30,000 participants.

## Run from one entry point

The repository includes the SNP list at `metadata/rsids.txt`, configured by the
`rsids` field in `config/pipeline.json`. On the **login node**, supply the ancestry
metadata, then run:

```bash
# Pilot: if panel preparation is needed, submit it and stop; otherwise submit LD.
bash scripts/run_pipeline.sh --chromosomes 22

# Production: reuse the prepared panel when its source list is unchanged.
bash scripts/run_pipeline.sh --chromosomes 1-22 --max-concurrent 4

# Validate existing dependencies and inputs without downloading or submitting.
bash scripts/run_pipeline.sh --chromosomes 22 --check-only

# Prepare references (may submit a panel job), without submitting LD.
bash scripts/run_pipeline.sh --chromosomes 22 --prepare-only
```

This entry point loads the Python/PLINK modules, sets up or checks the environment,
checks metadata and genotype headers, and performs small reference downloads on
the login node. **The full dbSNP scan, checksum verification, and rsID panel
construction run only in a SLURM job.** PVAR scans, full panel validation and LD
computation also run on compute nodes.

If the panel is missing or stale, the script ensures the VCF is local, submits
`scripts/prepare_panel.sbatch`, prints the job ID and a request to wait, and exits.
**No LD jobs are submitted in that invocation.** Wait for completion, check
`logs/panel-<job-id>.out`, then rerun the same command. For example:

```bash
bash scripts/run_pipeline.sh --chromosomes 22
# If a panel job was submitted, wait and inspect its progress:
squeue -j <job-id>
# After it finishes successfully:
bash scripts/run_pipeline.sh --chromosomes 22
```

Rerunning while that panel job is queued/running reports the existing job instead
of submitting it again, including when its source list has changed. If it fails,
inspect its log, fix the cause, and rerun; the script can resubmit preparation.
Panel-job state is stored under `logs/panel-<output-path-hash>.json`. If SLURM
status cannot be queried, the script fails rather than risk a duplicate job.

Default panel allocation: **1 CPU, 32 GiB RAM, 5 hours**, adjustable in
`prepare_panel.sbatch`. It uses the same Python module and environment as LD jobs.
All panel-job inputs must be on shared storage; it performs no downloads.
Only after a completed panel is available and the lightweight checks succeed
does a normal run submit the LD chromosome array. `--prepare-only` can submit a
panel preparation job but never LD jobs. `--check-only` never installs, downloads,
or submits anything; a missing/stale panel fails validation. Successful checks
write `logs/preflight.json` and check that log/output/configured scratch directories
are writable. Jobs otherwise rely on scheduler-provided `SLURM_TMPDIR`.

Use `--config FILE`, `--chromosomes 1-3,22`, and `--max-concurrent N` as needed.
`logs/preflight.json` records `awaiting_panel` or the validated LD submission/job ID.
LD jobs receive a snapshot of the validated configuration. `submit.sh` is an alias
for this workflow; direct `sbatch` bypasses login-node preflight.

The default rsID source is `metadata/rsids.txt`. Change `rsids` in the configuration
or use `--rsids FILE [FILE...]`. Login-node freshness checks compare file paths,
sizes and modification times against a completion receipt; they do not load the
rsID lists or hash the full panel. Input changes trigger another preparation job.
Hashes and per-chromosome variant counts are recorded by the compute job. Files
modified without changing size/mtime cannot be detected by these lightweight checks;
remove the panel receipt to force regeneration. Older receipts without panel-file
metadata require a one-time rebuild. Set `rsids` to `[]` only for an independently
prepared panel; its full validation is deferred to the chromosome jobs.

The default dbSNP location is `data/reference/GCF_000001405.40.gz`. If absent, the
login script downloads it with resume support and saves the accompanying MD5.
**The checksum scan is deferred to the panel job.** Use `--dbsnp-vcf FILE` to reuse
a local GRCh38 VCF, or `--dbsnp-url URL` for another release with `URL.md5`.
If a local `<VCF>.md5` exists, the compute job verifies it; otherwise independently
provided VCFs are assumed to have been verified upstream. A completed panel can
be reused without keeping the full source VCF.

Genotype checks inspect readability, PGEN magic bytes, PVAR headers/first records,
PSAM IDs, and ancestry intersections. They do **not** scan the full PGEN/PVAR or
prove all files' sample/variant counts match. PLINK and the chromosome processing
stage perform the full checks on compute nodes. A header-only preflight keeps the
login-node work bounded for the 455-million-variant source dataset.

The individual setup/preparation commands below remain available for debugging
and manual control. The unified entry point honors configured panel, map and block
paths; block files are prepared from the configured ancestry mapping and checked
against their provenance manifest before submission.

## Design and filtering plan

1. Supply a union of HapMap3 and your chosen PRS panels as one or more rsID lists.
   There is no universal modern PRS panel: HapMap3 alone is smaller than the desired
   2–2.5 million sites. Panel membership, GRCh38 availability, and ancestry QC
   determine the final size; do not LD-prune or randomly thin a reference matrix.
2. Stream a GRCh38 dbSNP VCF once to resolve **only requested rsIDs** into a compact
   allele-aware panel. Match chromosome, position, REF, and ALT exactly against PVAR;
   never map by position alone, strand-flip, or use GRCh37 coordinates. Keep only
   biallelic variants with distinct, nonempty A/C/G/T sequence alleles, including
   insertions/deletions such as `chr10:10370:A:AC` and multi-base substitutions.
   Symbolic/missing alleles and multiallelic records remain excluded. There is no
   arbitrary sequence-length cutoff. Reject ambiguous rsID/locus mappings and
   duplicate source records. rsID aliases absent from the supplied dbSNP
   release are reported as unresolved; reconcile merged IDs upstream if needed.
   Matching remains exact: differently padded/left-aligned indels are not
   automatically normalized or matched by position alone. Both sources must use
   the same GRCh38 representation. No reference genome is downloaded for normalization.
   Panel receipts record `variant_policy`; old SNP-only panels automatically trigger
   a SLURM rebuild, after which the launcher exits and asks you to wait and rerun.
   Independently supplied panels (`rsids: []`) remain under the user's control.
   The new summary distinguishes `rsid_outcomes` (disjoint unique-ID counts),
   `excluded_records` (VCF-row counts), and `resolved_by_type` (SNP/indel/MNV).
3. Validate ancestry metadata and intersect FID/IID pairs with each chromosome's
   PSAM. Write ancestry keep files and report unmatched labels. Optional sample QC
   keep files can restrict to an independently selected unrelated, high-quality
   sample set. The pipeline does not infer ancestry or remove relatives itself.
4. Scan each chromosome's PVAR, select the small target set, then run PLINK2
   `--extract` before BED conversion. This avoids materializing all 455 million
   variants as BED or repeatedly processing the full PGEN for every ancestry.
   `.pvar`, `.pvar.gz`, and `.pvar.zst` are supported (zstd executable for the latter).
5. Convert retained dosages to hardcalls with distance threshold 0.1. BED cannot
   preserve dosages. Re-read BED before missingness/frequency QC so those statistics
   describe exactly the hardcalls used for LD. This threshold is dosage distance,
   **not** a posterior genotype probability threshold.
6. Within each ancestry: sample missingness <= 0.05 on the panel, followed by variant
   missingness <= 0.05, MAF >= 0.001, MAC >= 20, and HWE p >= 1e-6 with `midp
   keep-fewhet`. `keep-fewhet` targets excess heterozygosity and reduces rejection
   from residual structure; review diagnostics for admixed groups. Do not apply
   pooled HWE or pooled MAF, which could lose variants common in minority groups.
   MAC 20 gives an effective minimum MAF near 20/(2N), often stricter than 0.1%
   for small groups. Boundaries are inclusive. Defaults are starting points to
   assess on chr22, not a guarantee of LD precision or panel size.
   MAC 10 can be evaluated by setting `mac` to `10`. Pearson correlations remain
   computable, but ten minor copies may be carried by only 5–10 people; rare–rare
   correlations can be sensitive to individual carriers and genotype errors.
   This is an exploratory retention/precision trade-off, not a universal accuracy
   guarantee. MAC remains 20 by default; MAF >= 0.001 still applies independently.
7. Optional imputation R2 filtering (suggested starting threshold 0.8) requires an
   explicitly named scalar PVAR INFO field. Missing/nonfinite values fail when
   enabled. No field name or quality availability is assumed; by default no
   imputation-quality filter is applied. Verify upstream quality before production.
8. Interpolate cM from the Beagle-distributed HapMap **GRCh38** map into BIM column
   3; use the nearest endpoint beyond map coverage and report those counts.
9. Compute ancestries sequentially in one job per chromosome, using magenpy's
   native BED backend (no PLINK1 required). Default: **block-diagonal LD, int8**,
   Zstandard, consolidated Zarr v2 ZIP. GRCh38 pyrho blocks are matched to each
   ancestry using the HGDP+1KG mappings below. Block truncation and int8
   quantization may affect numerical properties; downstream methods may need
   regularization. The interpolated common cM map remains in BIM metadata.

[PLINK2 filtering](https://www.cog-genomics.org/plink/2.0/filter),
[dosage hardcalls](https://www.cog-genomics.org/plink/2.0/input#dosage_import), and
[magenpy LD](https://shz9.github.io/magenpy/tutorials/compute_ld_matrices/) document
these interfaces. Per-ancestry SNP sets may differ. Preserve BIM alleles and SNP
order when matching summary statistics; retain palindromic SNPs and harmonize them
with alleles/frequencies downstream.

## Inputs and setup

Run commands from the repository root on the cluster. No genotype downloads are
performed. Default input prefix:

```text
$HOME/links/projects/def-sgravel/cartagene/flagship_array_imputed_w_cag_and_topmed_r3/PGEN/chr{chrom}.imputed_metaminimac_w_CaG_TOPMed_r3
```

Provide `data/metadata/ancestry_labels.tsv` (tab-separated; IDs read as strings):

```text
FID	IID	AFR	AMR	CSA	EAS	EUR	MID	Ancestry
1111	1111	0.0	0.15	0.0	0.0	0.85	0.0	EUR
```

The assigned `Ancestry` column is authoritative; proportions are not reclassified.
`OTH` samples (unknown ancestry) are excluded before sample selection. Other
unexpected labels and duplicate IDs fail validation. Missing chromosome samples are
reported; absent groups are skipped. Groups below 100 post-QC samples are skipped;
groups below 1,000 merit particular caution. For unrelated LD references, supply a
headered `FID IID` keep file through `sample_keep` in the config.

Create a persistent environment once **on the login node**:

```bash
bash scripts/setup_env.sh
```

Environment setup, reference downloads, panel preparation, and LD jobs each load
`python/3.13` before invoking Python. Setup uses that module's `python3` to create
`env/ld`; subsequent steps use `env/ld/bin/python` with the same module loaded.
To use another supported module (Python 3.11–3.13), export `PYTHON_MODULE` before
setup, preparation, and submission, for example `export PYTHON_MODULE=python/3.12`.
Recreate the environment when changing Python versions. Load any additional
compiler modules required by your cluster before setup.

Setup installs binary wheels for magenpy 0.2.5 and Zarr v2 in `env/ld`; it fails
if suitable wheels are unavailable rather than compiling on the login node. Jobs use `module load
plink` and resolve `plink2` or `plink` from that module, verifying it is PLINK v2.
No PLINK download is performed. This venv route avoids assuming an available Docker
image/tag or compute-node network access. It must be compatible with compute nodes.

## Reference downloads and scheduled preparation

Compute nodes have no internet access. Download public reference files on the
login node; keep them and `env/ld` on storage visible to compute nodes. Small
block/map tables can be normalized there. Full database scans and large panel
operations must run under SLURM. The unified entry point handles this split.

The bundled `metadata/rsids.txt` uses one rsID per line, no header (blank lines
and # comments allowed). For an existing GRCh38 dbSNP copy, use:

```bash
bash scripts/run_pipeline.sh --chromosomes 22 --dbsnp-vcf /shared/GRCh38.vcf.gz
```

The legacy preparation wrapper delegates to the unified entry point with
`--prepare-only`; it also schedules panel preparation and returns rather than
scanning dbSNP on login:

```bash
bash scripts/prepare_references.sh /shared/GRCh38.vcf.gz metadata/rsids.txt
```

For manual submission when all source files and the Python environment are local:

```bash
mkdir -p logs
sbatch scripts/prepare_panel.sbatch /shared/GRCh38.vcf.gz \
  data/reference/panel.grch38.tsv metadata/rsids.txt
```

`prepare_panel.py` itself refuses to run outside a SLURM job. Its JSON completion
receipt is published only after the panel succeeds. A failed or interrupted job
leaves no current receipt, so it cannot enable LD submission. Do not run heavy
checksum or panel validation commands on the login node.

LD blocks come from [LDblocks_GRCh38](https://github.com/jmacdon/LDblocks_GRCh38),
using **pyrho** boundaries. The `block_ancestry_map` object in
`config/pipeline.json` defines the source ancestry for every CARTaGENE ancestry;
the defaults match HGDP+1KG:

| CARTaGENE ancestry | Source block ancestry | Local normalized file |
|---|---|---|
| AFR | AFR | `data/blocks/AFR.bed` |
| AMR | AFR (proxy) | `data/blocks/AMR.bed` |
| CSA | SAS (proxy) | `data/blocks/CSA.bed` |
| EAS | EAS | `data/blocks/EAS.bed` |
| EUR | EUR | `data/blocks/EUR.bed` |
| MID | EUR (proxy) | `data/blocks/MID.bed` |

`prepare_blocks.py` reads these mappings from `pipeline.json`; they are not
hardcoded in Python. Edit the mapping and rerun login-node reference preparation
to regenerate affected files. Outputs retain the CARTaGENE ancestry names, so
`blocks: "data/blocks/{ancestry}.bed"` continues to work for chromosome jobs.
For a different configuration, use `prepare_blocks.py --config /path/to/pipeline.json`
or `bash scripts/prepare_references.sh --config /path/to/pipeline.json VCF panel...`.
The preparation wrapper honors the unified workflow’s configured reference paths.

The proxies preserve the previous pipeline's choices; they are not independently
estimated CARTaGENE or admixed-population boundaries. Source coordinates are
preserved and passed to magenpy as `[start, stop)` intervals. SNPs in gaps or
outside the blocks are excluded from LD; QC and LD SNP counts are reported
separately. Cross-block correlations are set to zero.

Block preparation validates non-overlapping intervals and coverage of all 22
autosomes. `data/blocks/summary.json` records each ancestry's source, URL, counts,
and output SHA256. Cached blocks are reused only when their manifest matches.
Download/validation errors stop preparation rather than silently dropping groups.
The map archive and normalized maps are cached in `data/genetic_map/`, with source
metadata in `summary.json`. Panel preparation records dbSNP file size/mtime and
panel-list hashes; `--assembly GRCh38` is a required assertion. To refresh blocks,
run `env/ld/bin/python scripts/prepare_blocks.py --force` on the login node.

To opt into windowed LD, set `estimator` to `windowed` and `window_cm` to `3` in
the configuration. The common cM map is an approximation for all groups.
`storage_dtype` can independently be changed to `int16` or floating point.

## Pilot and production

Edit `config/pipeline.json` for paths, groups, QC and output settings. Relative
paths resolve from the repository root; environment variables and `~` expand.
The unified entry point snapshots the validated configuration for submitted jobs.
Keep reference and genotype files fixed while an array runs.

```bash
bash scripts/run_pipeline.sh --config config/pipeline.json --chromosomes 22
# Inspect counts, PLINK logs, map endpoint counts and LD before scaling:
bash scripts/run_pipeline.sh --config config/pipeline.json --chromosomes 1-22 --max-concurrent 4
```

Default allocation per chromosome: 16 CPUs, 128 GiB RAM, 48 hours; at most four
chromosomes concurrently in the example. These are initial requests, not measured
requirements. chr1 is a second useful resource pilot. Use `sacct` to inspect peak
RSS and runtime and adjust `--mem`/`--time` via `SBATCH_MEM`/`SBATCH_TIMELIMIT` or
edit the batch directives. Large blocks/high-density regions can dominate storage
and compute. Shared full-panel BED is about N*M/4 bytes (roughly 19 GB genome-wide
at 30k x 2.5M), excluding per-ancestry copies, LD and temporaries. Provision generous
scratch; LD size depends on local density and block size (or map plateaus in window mode). Jobs require
`SLURM_TMPDIR` or a configured `scratch_dir` on a suitable filesystem.

## Outputs and reruns

Each job fingerprints configuration, scripts, reference/metadata content, input
PGEN/PVAR/PSAM size+mtime, and tool versions. Completed results live under:

```text
results/<chromosome-fingerprint>/chr22/
  summary.json
  commands.jsonl
  shared.*.log
  EUR/
    samples.keep
    qc.bim, qc.fam, qc.log, samples.log
    variants.tsv
    ld/chr_22.zip
```

The exact output path is printed in the job log. Summary reports panel matches,
unmatched metadata, pre/post-QC sample/variant counts, LD variant counts, skipped groups, endpoint cM
counts, parameters and software versions. `variants.tsv` records original ID,
rsID, GRCh38 REF/ALT and position. Keep BIM/FAM and logs private: they contain sample
IDs. Data, environments, logs and results are gitignored. No upload or publication
is performed. QC BEDs are copied only when `retain_bed` is true; scratch is cleaned
on normal exit, and the scheduler manages node-local scratch after termination.

A success receipt is written only after ZIP CRC and consolidated-metadata checks
and atomic publication of a chromosome directory. Existing receipts with intact
archives skip recomputation. Incomplete jobs are recomputed in fresh scratch;
outputs for changed inputs use a new fingerprint. Large PGENs are identified by
size/mtime rather than fully hashed; replace files with updated mtimes or change
`run_tag` when their content changes. Do not submit duplicate chromosome tasks to
the same output configuration concurrently. An output lock prevents collisions;
a killed job may leave a lock requiring manual removal after confirming it is dead.

## Local verification

```bash
PYTHONDONTWRITEBYTECODE=1 env/ld/bin/python -m unittest discover -s tests -v
for script in scripts/*.sh scripts/*.sbatch; do bash -n "$script"; done
```

Pandas and NumPy (installed with magenpy) handle metadata tables and interpolation;
full dbSNP and PVAR files remain streamed to bound memory. Tests cover
allele-aware mapping, ambiguous sites, sample intersections and cM
annotation with synthetic data. When magenpy is installed, an additional numerical
test checks windowed and ancestry-specific block correlations against NumPy,
including int8 storage, cross-block zeros and reopening saved ZIPs; otherwise
that test is skipped. Run the full set with `env/ld/bin/python -m unittest discover
-s tests -v`. To inspect a saved matrix with magenpy:

```python
import magenpy as mgp
ld = mgp.LDMatrix.from_path("results/<fingerprint>/chr22/EUR/ld/chr_22.zip")
print(ld.shape, ld.sample_size)
# For a SMALL matrix only: float dtype dequantizes; symmetry restores the diagonal.
r = ld.to_csr(dtype="float32", return_symmetric=True)
```

End-to-end computation requires cluster data,
PLINK2 and the configured magenpy environment; run the chr22 pilot before scaling.
