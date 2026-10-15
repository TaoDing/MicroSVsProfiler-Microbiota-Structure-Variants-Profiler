# Step 2 - Microbial taxonomic and functional profiling

## Purpose

This step profiles microbial taxonomy with **Kraken2** and **Bracken** and functional potential with **HUMAnN 3.0**. Inputs are quality-controlled, host-depleted paired-end FASTQ files from the preceding preprocessing step.

**Terminology:** the workflow yields taxonomic abundance estimates and gene family/pathway profiles. These outputs alone do not establish microbial gene expression or measured metabolic activity.

## 1. Software and databases

| Component | Version reported in the original procedure | Purpose |
| --- | --- | --- |
| Kraken2 | 2.1.2 | Read classification |
| Bracken | 2.6.2 | Abundance re-estimation, species level (`-l S`) |
| HUMAnN | 3.0 | Gene family and pathway profiling |
| UniRef90 | UniProt release 2019_01 | Protein reference database for HUMAnN |
| MetaPhlAn database | mpa_v30_CHOCOPhlAn_201901 | Taxonomic prescreening in HUMAnN |
| Kraken2 database | 5/17/2021 | Kraken2/Bracken classification |

Confirm installed package versions and database releases from the computing environment. In particular, the exact Kraken2 database build, MetaPhlAn database/index, host reference version (previous step), and HUMAnN database releases are important for reproducibility.

Verify installations:

```bash
kraken2 --version
bracken -v  # if unsupported in your installation, use bracken -h
humann --version
metaphlan --version
bowtie2 --version
```

Bracken requires a database-specific k-mer/read-length distribution generated for the requested `-r` setting. The procedure uses `-r 150`; select the value that matches the prepared Bracken distribution. Trimmed reads may have variable lengths; document the chosen configuration.

## 2. Inputs

The same biological samples are used by both analysis branches. Create a tab-separated sample manifest with these **exact** column names:

```tsv
sample_id	read1	read2
SAMPLE001	/path/to/step1/SAMPLE001_R1.fastq.gz	/path/to/step1/SAMPLE001_R2.fastq.gz
```

See [`examples/profiling_samples.example.tsv`](../examples/profiling_samples.example.tsv). Each row must identify a pair of synchronized FASTQ files. Uncompressed `.fastq` and gzip-compressed `.fastq.gz` are supported; mates must have matching compression. The files come from the actual Step 1 host-read removal and subsampling outputs, and the manifest paths must be edited accordingly.

Specify directories for:

- a complete **Kraken2 database** compatible with the Bracken database files;
- **MetaPhlAn Bowtie2 database** matching the installed HUMAnN/MetaPhlAn release;
- HUMAnN **ChocoPhlAn nucleotide** database;
- HUMAnN **UniRef90 protein** database.

These databases are not included in this repository. Their provenance and precise release versions should be recorded.

## 3. Execute taxonomic and functional profiling

From the repository root:

```bash
bash scripts/step2_taxonomic_functional_profiling.sh \
  --samples examples/profiling_samples.example.tsv \
  --kraken-db /path/to/kraken2_database \
  --metaphlan-db /path/to/metaphlan_database \
  --humann-nucleotide-db /path/to/chocophlan_database \
  --humann-protein-db /path/to/uniref90_database \
  --out results/step2 \
  --threads 20 \
  --bracken-read-length 150
```

The script executes each sample sequentially, allocating the requested number of threads to tools. For cluster execution, submit the command as a scheduler job after adapting CPU/memory requests to local policies; it does not assume a PBS queue name. Run one manifest/subset per independent parallel job if necessary. Available disk space should accommodate HUMAnN intermediate FASTQ files and alignment results.

### 3.1 Kraken2 classification

Kraken2 classifies quality-controlled paired reads using the supplied database and emits a per-read classification and taxonomic report:

```bash
kraken2 --db /path/to/kraken2_database --threads 20 \
  --paired --gzip-compressed \
  --report results/step2/kraken2/SAMPLE001.report \
  --output results/step2/kraken2/SAMPLE001.kraken.output \
  SAMPLE001_R1.fastq.gz SAMPLE001_R2.fastq.gz
```

The executable script automatically omits `--gzip-compressed` for uncompressed input.

### 3.2 Bracken species-abundance estimation

```bash
bracken -d /path/to/kraken2_database \
  -i results/step2/kraken2/SAMPLE001.report \
  -o results/step2/bracken/SAMPLE001.bracken.tsv \
  -w results/step2/bracken/SAMPLE001.bracken.report \
  -r 150 -l S
```

Bracken reports estimated species-level read counts and associated abundance fractions. For downstream analyses, a sample-by-species matrix can be constructed from the Bracken reports, with relative abundance calculated consistently for the chosen set of species. The exact denominator and feature filtering should be documented in the downstream analysis method; the current script does **not** silently apply extra filters.

### 3.3 Prepare HUMAnN input

```bash
python3 scripts/merge_paired_fastq.py \
  --r1 SAMPLE001_R1.fastq.gz \
  --r2 SAMPLE001_R2.fastq.gz \
  --output results/step2/humann_input/SAMPLE001.fastq
```

The Python program validates FASTQ record structure, matching read identifiers, and equal mate counts. It writes both mate reads into a single FASTQ file in consecutive R1/R2 order. **This is read concatenation for HUMAnN input, not overlap-based read merging or a reconstruction of longer DNA molecules.**

### 3.4 HUMAnN functional profiling

```bash
humann --input results/step2/humann_input/SAMPLE001.fastq \
  --output results/step2/humann/SAMPLE001 \
  --threads 20 \
  --nucleotide-database /path/to/chocophlan_database \
  --protein-database /path/to/uniref90_database \
  --metaphlan-options "--bowtie2db /path/to/metaphlan_database"
```

HUMAnN profiles gene families and metabolic pathways with its database-dependent alignment workflow. Its outputs typically include `*_genefamilies.tsv`, `*_pathabundance.tsv`, and `*_pathcoverage.tsv` (alongside logs and intermediate files). Read-equivalent/RPK-style units and normalization depend on HUMAnN output stage and chosen normalization commands; **this script does not automatically normalize all outputs to CPM or relative abundance**. If required, use `humann_renorm_table` and document the selected `--units` option. The exact MetaPhlAn release/index must be compatible with the deployed HUMAnN environment.

## 4. Expected output tree

```text
results/step2/
├── kraken2/
│   ├── SAMPLE001.report
│   └── SAMPLE001.kraken.output
├── bracken/
│   ├── SAMPLE001.bracken.tsv
│   └── SAMPLE001.bracken.report
├── humann_input/
│   └── SAMPLE001.fastq
├── humann/
│   └── SAMPLE001/
│       ├── SAMPLE001_genefamilies.tsv
│       ├── SAMPLE001_pathabundance.tsv
│       └── SAMPLE001_pathcoverage.tsv
└── logs/
    ├── SAMPLE001.kraken.stderr.log
    ├── SAMPLE001.bracken.log
    ├── SAMPLE001.merge.log
    └── SAMPLE001.humann.stderr.log
```

Output filenames for HUMAnN can vary with input basename and HUMAnN version. Inspect the actual output directory rather than assuming names without checking.

## 5. Quality-control checks and reproducibility notes

1. Confirm each manifest ID maps uniquely to one synchronized paired-end sample, with no missing mates.
2. Examine Kraken2 classification summaries, species reports, and Bracken reports for plausible read totals.
3. Check the Bracken `-r` database distribution and the Kraken2 database build/release.
4. Check HUMAnN log files, mapped/unmapped read fractions, and availability/compatibility of ChocoPhlAn, UniRef90, and MetaPhlAn resources.
5. Record the exact software versions, reference database releases, read preprocessing and subsampling parameters, and the reference database paths used on the analysis server.
6. Document the method used to combine per-sample tables and normalize species/pathway abundances in subsequent analyses.