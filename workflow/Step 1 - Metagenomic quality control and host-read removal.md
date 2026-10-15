# Metagenomic quality control and host-read removal

This workflow document described the two-stage workflow in the analysis: (1) quality control and host-read removal with KneadData, and (2) paired-end read subsampling with seqtk. Execution is implemented in Bash and Python.

## 1. Software and reference requirements

| Component | Version | Purpose |
| --- | --- | --- |
| KneadData | 0.7.5 | Trimming and host-read filtering |
| Trimmomatic | 0.39 | Quality trimming |
| Bowtie2 | 2.2.3 | Mapping to host reference |
| seqtk | 1.3 | Random read-pair subsampling |

## 2. Input files

Prepare a tab-delimited sample manifest such as [`examples/samples.example.tsv`](../examples/samples.example.tsv):

```tsv
sample_id	read1	read2
SAMPLE001	/path/to/SAMPLE001_1.fastq.gz	/path/to/SAMPLE001_2.fastq.gz
```

- Exactly one sample per line; first line contains column names.
- Both mate files must contain read pairs with matching identifiers and order.
- FASTQ files may be plain `.fastq` or gzip-compressed `.fastq.gz`.
- `read1` and `read2` may be absolute paths or paths relative to the manifest directory.
- Provide a prepared Bowtie2-compatible human reference **directory** as the second script argument.

## 3. Run

Install and activate the environment containing `kneaddata`, `seqtk`, and `python3`. From the repository root:

```bash
bash scripts/step1_preprocess.sh \
  examples/samples.example.tsv \
  /path/to/kneaddata_human_reference/ \
  results/step1 \
  20
```

Parameters are `manifest`, `host_reference_dir`, `output_dir`, and optional `threads` (default 20).

### Running on a PBS cluster

Run one sample per PBS job by preparing a one-sample manifest, then executing the command above within the job script. Example resource directives (change for your scheduler):

```bash
#!/bin/bash
#PBS -N MicroSV_Step1
#PBS -l nodes=1:ppn=20
#PBS -j oe
set -euo pipefail
cd /path/to/MicroSVsProfiler
# Activate the appropriate conda environment here.
bash scripts/step1_preprocess.sh \
  /path/to/one_sample.tsv /path/to/human_reference \
  results/step1 20
```

Submit using `qsub job_step1.pbs` where supported. The supplied Bash script runs samples sequentially; it does not submit PBS jobs automatically.

## 4. Outputs

```text
results/step1/
├── kneaddata/
│   └── SAMPLE001/
│       ├── ..._kneaddata_paired_1.fastq
│       ├── ..._kneaddata_paired_2.fastq
│       └── SAMPLE001.kneaddata.log
└── seqtk/
    ├── SAMPLE001_R1.fastq
    └── SAMPLE001_R2.fastq
```

KneadData will produce additional intermediate and unpaired files; only its paired, filtered reads are downsampled. Exact KneadData paired-output names depend on the input naming convention and installed version. The wrapper explicitly detects one paired file per mate and fails if the expected outputs are missing or ambiguous. The FASTQ validator rejects unmatched read identifiers, differing mate counts, and malformed four-line records.

## 5. Reproducibility checklist

Record the exact tool and host-reference versions, Trimmomatic options actually executed, database download/accession, sample manifest, random seed, and command-line logs. Consider validating the pipeline on a small sample before starting a full-cohort run.
