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
# Pilot: prepare/check dependencies and inputs, then submit chromosome 22.
bash scripts/run_pipeline.sh --chromosomes 22

# Production: reuse the prepared panel when its source list is unchanged.
bash scripts/run_pipeline.sh --chromosomes 1-22 --max-concurrent 4

# Validate existing dependencies and inputs without downloading or submitting.
bash scripts/run_pipeline.sh --chromosomes 22 --check-only

# Prepare and validate everything, stopping before submission.
bash scripts/run_pipeline.sh --chromosomes 22 --prepare-only
```

This entry point loads the Python/PLINK modules, creates or repairs the Python
environment when needed, checks installed dependencies and the compiled native
LD backend, validates metadata and **all requested genotype triplets**, prepares
missing references, validates them, and submits a SLURM array only after every
check succeeds. A failure prevents submission of the entire requested array.
`--check-only` still writes a validation report to `logs/preflight.json` and checks
that output/log/scratch directories are writable; it installs/downloads nothing.
An explicitly configured scratch directory is checked on the login node;
otherwise jobs rely on the scheduler-provided `SLURM_TMPDIR`.

Use `--config FILE` for alternate settings, `--chromosomes 1-3,22` to select
chromosomes, and `--max-concurrent N` to limit simultaneous jobs. Normal successful
runs **submit jobs automatically**. `logs/preflight.json` records validation details
and the submitted job ID. Jobs receive a snapshot of the validated configuration.
`submit.sh` remains an alias for this workflow; direct `sbatch` is an advanced
path that bypasses the login-node checks.

The default rsID source is `metadata/rsids.txt`; no `--rsids` argument is needed.
Change the `rsids` list in the configuration or pass `--rsids FILE [FILE...]` to
override it. The source list is checked against the panel receipt on each run;
changes trigger a panel rebuild (or a stale-panel error in `--check-only` mode).
Set `rsids` to `[]` only when supplying an independently prepared panel without
source lists. The default dbSNP location is
`data/reference/GCF_000001405.40.gz`. If absent, the script downloads it on the
login node with resume support and verifies the source's accompanying MD5 before
use. Set `--dbsnp-vcf /path/to/existing/GRCh38.vcf.gz` to reuse a verified local copy,
or `--dbsnp-url URL` to choose a GRCh38 release (with an accompanying `URL.md5`).
Existing local VCFs are assumed to have been verified by you. A supplied panel is
validated without requiring the full VCF; when source lists are configured or supplied, their
hashes, the panel receipt, and available VCF metadata determine whether to rebuild.

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
   biallelic A/C/G/T SNPs. Reject ambiguous rsID/locus mappings and duplicate source
   records. Indels such as `chr10:10370:A:AC` are intentionally excluded in this
   initial SNP-panel implementation. rsID aliases absent from the supplied dbSNP
   release are reported as unresolved; reconcile merged IDs upstream if needed.
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
6. Within each ancestry: sample missingness <= 0.02 on the panel, followed by variant
   missingness <= 0.02, MAF >= 0.001, MAC >= 20, and HWE p >= 1e-6 with `midp
   keep-fewhet`. `keep-fewhet` targets excess heterozygosity and reduces rejection
   from residual structure; review diagnostics for admixed groups. Do not apply
   pooled HWE or pooled MAF, which could lose variants common in minority groups.
   MAC 20 gives an effective minimum MAF near 20/(2N), often stricter than 0.1%
   for small groups. Boundaries are inclusive. Defaults are starting points to
   assess on chr22, not a guarantee of LD precision or panel size.
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

Environment setup, reference preparation, and chromosome jobs each load
`python/3.13` before invoking Python. Setup uses that module's `python3` to create
`env/ld`; subsequent steps use `env/ld/bin/python` with the same module loaded.
To use another supported module (Python 3.11–3.13), export `PYTHON_MODULE` before
setup, preparation, and submission, for example `export PYTHON_MODULE=python/3.12`.
Recreate the environment when changing Python versions. Load any additional
compiler modules required by your cluster before setup.

This installs pinned magenpy 0.2.5 and Zarr v2 in `env/ld`. Jobs use `module load
plink` and resolve `plink2` or `plink` from that module, verifying it is PLINK v2.
No PLINK download is performed. This venv route avoids assuming an available Docker
image/tag or compute-node network access. It must be compatible with compute nodes.

## Prepare public references on the login node

**Compute nodes have no internet access.** Run environment setup, all downloads,
and reference preparation directly on the login node before submitting LD jobs.
Do not use `sbatch` for these steps. The preparation entry points reject invocation
inside a SLURM job. There is no download fallback in the chromosome jobs: missing
local references cause an early failure with preparation instructions.

The bundled `metadata/rsids.txt` contains one rsID per line, without a header
(blank lines and # comments are allowed). Additional lists can be selected through
`rsids` in the configuration or `--rsids`. Obtain a GRCh38
dbSNP VCF from [NCBI](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/) on the
login node, record its release/checksum, and verify the assembly. This is a large
download; reuse the verified HGDP+1KG copy if available. For example, on the login node:

```bash
mkdir -p data/reference
curl --fail --location --continue-at - \
  --output data/reference/GCF_000001405.40.gz \
  https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
# Verify the downloaded release against NCBI's accompanying checksum before use.
bash scripts/prepare_references.sh data/reference/GCF_000001405.40.gz \
  metadata/rsids.txt
```

The login-node script downloads and normalizes the Beagle GRCh38 genetic map and
ancestry-matched LD blocks, then streams the local dbSNP VCF to build
`data/reference/panel.grch38.tsv`. The full dbSNP scan may take hours; allocate
sufficient login-node resources according to cluster policy. This replaces the
former `prepare_references.sbatch` workflow. All reference outputs and `env/ld`
must be on storage visible to compute nodes.

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
The preparation wrapper's other reference destinations remain the standard paths.

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
