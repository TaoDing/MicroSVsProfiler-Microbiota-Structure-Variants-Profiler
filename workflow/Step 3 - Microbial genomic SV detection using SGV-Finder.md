# Step 3 - Microbial genomic SV detection using SGV-Finder

## Overview

This stage detects **coverage-defined microbial genomic structural variation (SV) features** from quality-controlled shotgun metagenomic reads. It follows the general strategy of Zeevi et al. (2019) and the original SGV-Finder pipeline:

1. Map reads to representative microbial genomes using GEM Mapper and refine ambiguous assignments with **ICRA**.
2. Calculate per-sample genomic coverage using **SGVF_PerFile**.
3. Aggregate samples **by cohort**, identify deletion-like SV regions (**dSVs**) and highly variable coverage regions (**vSVs**) using **SGV-Finder**, and write results to CSV.

**Important interpretation:** These are coverage-defined genomic regions, *not* experimentally resolved insertions, deletions, inversions or breakpoint coordinates. SV feature IDs should not be interpreted as physical lengths or positions.

## Prerequisites

- Step 1 output: quality-controlled, host-depleted paired-end reads (and, if applicable, subsampled reads).
- SGV-Finder/ICRA source tree, compatible Python environment and **GEM Mapper** installed according to the software's original instructions.
- SGV-Finder reference data: 3,953 representative proGenomes-derived microbial genomes, as stated in the original protocol. Record the exact source, archive checksum and local path when running the analysis.
- Sufficient memory and compute resources for read mapping and cohort-wide coverage processing.

> **Legacy compatibility:** The original protocol activated a Python 2.7.9 environment for the SGV-Finder vendor commands. This wrapper uses Bash, and the optional export script uses Python 3, but **does not rewrite SGV-Finder itself**. Set `--python` to the interpreter supported by your installed revision; do not assume current Python 3 works with the legacy vendor code.

## Files and directory layout

From the repository root:

```text
MicroSVsProfiler/
├── workflow/
│   └── Step 3 - Microbial genomic SV detection using SGV-Finder.md
├── scripts/
│   ├── step3_sgv_pipeline.sh
│   └── convert_sgv_pickle_to_csv.py
├── examples/
│   ├── step3_samples.example.tsv
│   └── step3_cohorts.example.tsv
└── results/
    └── step3/
        ├── icra/<sample_id>/
        ├── perfile/<sample_id>.PerFile.jsdel
        ├── sv/<cohort>/<cohort>.dsgv.csv
        ├── sv/<cohort>/<cohort>.vsgv.csv
        └── converted/{dsgv.csv,vsgv.csv}  # optional
```

The two example TSV files are *templates*, not executable manifests. Replace placeholders with actual paths after checking the outputs and formats of your installed programs.

## A. Per-sample ICRA read reassignment

The ICRA command's input is a **read-file prefix** understood by the installed ICRA implementation, not the first FASTQ path alone.

```bash
bash scripts/step3_sgv_pipeline.sh icra \
  --sgv-src /path/to/SGVFinder/src \
  --sample SAMPLE_A \
  --reads-prefix /path/to/clean/SAMPLE_A_1_kneaddata_paired \
  --python /path/to/compatible/python \
  --out results/step3
```

**Verify:** The ICRA call and its path conventions must be checked against the installed source before execution. After successful execution, locate the actual generated `.jsdel` file and record its path; do not guess the output filename.

## B. Per-sample SGVF_PerFile coverage profiling

Used a fixed effective read length of **150**, `--x_coverage 0.01`, and `--rate_param 10`:

```bash
bash scripts/step3_sgv_pipeline.sh perfile \
  --sgv-src /path/to/SGVFinder/src \
  --sample SAMPLE_A \
  --jsdel /path/to/actual/SAMPLE_A_ICRA_output.jsdel \
  --read-length 150 \
  --python /path/to/compatible/python \
  --out results/step3
```

For reads of other lengths, verify the correct parameter with the installed SGVF_PerFile documentation.

## C. Cohort-level SV calling

**Run the final SGV-Finder command separately for each cohort**, on *multiple per-sample* `PerFile.jsdel` outputs.

```bash
bash scripts/step3_sgv_pipeline.sh call \
  --sgv-src /path/to/SGVFinder/src \
  --cohort COHORT_A \
  --input-glob 'results/step3/perfile/COHORT_A_*.PerFile.jsdel' \
  --min-samp-cutoff 75 \
  --python /path/to/compatible/python \
  --out results/step3
```

**Do not** supply a glob containing samples from unrelated cohorts. Ensure each cohort contains at least the appropriate number of qualifying samples, and check whether `--min_samp_cutoff` is accepted by your exact SGV-Finder version. A species present in fewer than 75 adequately covered samples cannot satisfy this parameter. A `5%` prevalence claim in the historical document should not be confused with `min_samp_cutoff = 75`; the latter is an absolute sample threshold.

## D. Optional pandas-format export

If your SGV-Finder implementation instead produces trusted pandas `dsgv.df` and `vsgv.df` pickle files, use:

```bash
python3 scripts/convert_sgv_pickle_to_csv.py \
  --dsv /path/to/dsgv.df \
  --vsv /path/to/vsgv.df \
  --out results/step3/converted
```

**Security:** `pandas.read_pickle()` can execute arbitrary code embedded in malicious pickle files. Only convert files generated by your own trusted analysis.

## Parameters and methodological context

| Item | Description |
|---|---|
| Reference genomes | 3,953 representative microbial genomes from proGenomes |
| Binning | 1-kb genomic bins |
| Genome coverage screening | Median ≥10 mapped reads per bin |
| Cohort eligibility | `min_samp_cutoff = 75` |
| dSV definition | Deletion-like regions in 25–75% of adequately covered samples |
| vSV definition | Substantial between-sample variation in normalized coverage |
| PerFile read length | `150` |
| Coverage setting | `--x_coverage 0.01` |
| Rate setting | `--rate_param 10` |
| Output options | `--csv_output --byorig` |

## Validation and reproducibility checks

1. Run `bash -n scripts/step3_sgv_pipeline.sh` and `python3 -m py_compile scripts/convert_sgv_pickle_to_csv.py`.
2. Confirm the installed program names and `--help` options before attempting the pipeline.
3. Record exact software revisions (`git rev-parse HEAD`), reference database release, paths, and checksums.
4. Check per-sample ICRA and PerFile logs and confirm each expected sample appears exactly once in its cohort.
5. Check the generated dSV/vSV matrix orientation (features versus samples), feature-ID definitions, and number of eligible species before interpreting results.

## References

- Zeevi D, et al. Structural variation in the gut microbiome associates with host health. *Nature*. 2019. https://doi.org/10.1038/s41586-019-1065-y
- SGV-Finder source: https://github.com/segalab/SGVFinder
