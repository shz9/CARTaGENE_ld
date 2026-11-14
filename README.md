# CARTaGENE ancestry-stratified GRCh38 LD

This pipeline turns chromosome-level imputed PGEN files into ancestry-specific,
QC-filtered BED files and signed LD correlation matrices using magenpy. Genotypes
stay on the SLURM cluster. Autosomes 1–22 are supported. The implementation follows
[HGDP+1KG](https://github.com/shz9/1KG_HGDP_ld), with early panel extraction to
handle CARTaGENE's approximately 455 million input variants and 30,000 participants.

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
   native BED backend (no PLINK1 required). Default: 3 cM window, int16, Zstandard,
   consolidated Zarr v2 ZIP. The common map is a distance approximation for all
   groups. Optional GRCh38 ancestry-specific blocks can replace windows; supply
   validated files for every requested ancestry. Window truncation and quantization
   may affect numerical properties; downstream methods may need regularization.

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
Unknown labels and duplicate IDs fail validation. Missing chromosome samples are
reported; absent groups are skipped. Groups below 100 post-QC samples are skipped;
groups below 1,000 merit particular caution. For unrelated LD references, supply a
headered `FID IID` keep file through `sample_keep` in the config.

Create a persistent environment once on a node with internet access (Python 3.11–3.13;
load your cluster's Python/compiler modules first if needed):

```bash
bash scripts/setup_env.sh
```

This installs pinned magenpy 0.2.5 and Zarr v2 in `env/ld`. Jobs use `module load
plink` and resolve `plink2` or `plink` from that module, verifying it is PLINK v2.
No PLINK download is performed. This venv route avoids assuming an available Docker
image/tag or compute-node network access. It must be compatible with compute nodes.

## Prepare public references once

Put your panel lists at e.g. `data/panels/hm3.rsids` and `data/panels/prs.rsids`,
one rsID per line, no header (blank lines and # comments allowed). Obtain a GRCh38
dbSNP VCF from [NCBI](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/), record its
release/checksum, and download it on a transfer node. This is a large download;
reuse the verified HGDP+1KG copy if available. The reference job streams it once:

```bash
mkdir -p logs
sbatch scripts/prepare_references.sbatch /absolute/path/to/GRCh38.dbsnp.vcf.gz \
  data/panels/hm3.rsids data/panels/prs.rsids
```

The job builds `data/reference/panel.grch38.tsv` and its JSON report, and downloads
and normalizes the Beagle GRCh38 map. It requires internet only if the map ZIP is
not already cached; to work offline pre-stage `data/genetic_map/plink.GRCh38.map.zip`.
It records input file sizes/mtimes and list hashes; archive source metadata is in
`data/genetic_map/summary.json`. Verify the dbSNP assembly independently; filenames
alone do not establish build. `--assembly GRCh38` is a required assertion.

## Pilot and production

Edit `config/pipeline.json` for paths, groups, QC and output settings. Relative
paths resolve from the repository root; environment variables and `~` expand.
The configuration is read at job start; keep it fixed while an array runs.

```bash
bash scripts/submit.sh config/pipeline.json 22
# Inspect counts, PLINK logs, map endpoint counts and LD before scaling:
bash scripts/submit.sh config/pipeline.json 1-22%4
```

Default allocation per chromosome: 16 CPUs, 128 GiB RAM, 48 hours; at most four
chromosomes concurrently in the example. These are initial requests, not measured
requirements. chr1 is a second useful resource pilot. Use `sacct` to inspect peak
RSS and runtime and adjust `--mem`/`--time` via `SBATCH_MEM`/`SBATCH_TIMELIMIT` or
edit the batch directives. Large windows/high-density regions can dominate storage
and compute. Shared full-panel BED is about N*M/4 bytes (roughly 19 GB genome-wide
at 30k x 2.5M), excluding per-ancestry copies, LD and temporaries. Provision generous
scratch; LD size depends on local density and map plateaus. Jobs require
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
unmatched metadata, pre/post-QC sample/variant counts, skipped groups, endpoint cM
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
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
for script in scripts/*.sh scripts/*.sbatch; do bash -n "$script"; done
```

Tests cover allele-aware mapping, ambiguous sites, sample intersections and cM
annotation with synthetic data. When magenpy is installed, an additional numerical
test checks signed correlations against NumPy and reopens the saved ZIP; otherwise
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
