#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/step2_taxonomic_functional_profiling.sh \
    --samples examples/profiling_samples.example.tsv \
    --kraken-db /path/to/kraken2_db \
    --metaphlan-db /path/to/metaphlan_db \
    --humann-nucleotide-db /path/to/chocophlan \
    --humann-protein-db /path/to/uniref90 \
    --out results/step2 [--threads 20] [--bracken-read-length 150]

Manifest must be tab-delimited with header: sample_id<TAB>read1<TAB>read2
EOF
}

samples='' kraken_db='' metaphlan_db='' nucleotide_db='' protein_db='' out='' threads=20 read_length=150
while (($#)); do
  case "$1" in
    --samples) samples=$2; shift 2;;
    --kraken-db) kraken_db=$2; shift 2;;
    --metaphlan-db) metaphlan_db=$2; shift 2;;
    --humann-nucleotide-db) nucleotide_db=$2; shift 2;;
    --humann-protein-db) protein_db=$2; shift 2;;
    --out) out=$2; shift 2;;
    --threads) threads=$2; shift 2;;
    --bracken-read-length) read_length=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done
for v in samples kraken_db metaphlan_db nucleotide_db protein_db out; do
  [[ -n ${!v} ]] || { echo "Missing --${v//_/-}" >&2; exit 2; }
done
[[ -f "$samples" ]] || { echo "Missing manifest: $samples" >&2; exit 2; }
for v in kraken_db metaphlan_db nucleotide_db protein_db; do
  [[ -d ${!v} ]] || { echo "Missing directory $v: ${!v}" >&2; exit 2; }
done
for app in kraken2 bracken humann python3; do
  command -v "$app" >/dev/null || { echo "Missing command: $app" >&2; exit 2; }
done
[[ $threads =~ ^[1-9][0-9]*$ && $read_length =~ ^[1-9][0-9]*$ ]] || { echo 'threads and read length must be positive integers' >&2; exit 2; }
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out"/{kraken2,bracken,humann_input,humann,logs}
header=$(head -n 1 "$samples" | tr -d '\r')
[[ $header == $'sample_id\tread1\tread2' ]] || { echo 'Manifest header must be sample_id<TAB>read1<TAB>read2' >&2; exit 2; }
# A manifest with an additional blank line at the end is accepted.
while IFS=$'\t' read -r sample r1 r2 extra || [[ -n ${sample:-} ]]; do
  sample=${sample%$'\r'}; r1=${r1%$'\r'}; r2=${r2%$'\r'}
  [[ -n $sample ]] || continue
  [[ $sample =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Unsafe sample ID: $sample" >&2; exit 2; }
  [[ -f "$r1" && -f "$r2" ]] || { echo "Missing reads for $sample: $r1 $r2" >&2; exit 2; }
  [[ -z ${extra:-} ]] || { echo "Extra columns for $sample" >&2; exit 2; }
  echo "[$sample] Kraken2" >&2
  compression=()
  if [[ $r1 == *.gz && $r2 == *.gz ]]; then compression=(--gzip-compressed)
  elif [[ $r1 == *.gz || $r2 == *.gz ]]; then echo 'Both mates must use matching compression' >&2; exit 2
  fi
  kraken2 --db "$kraken_db" --threads "$threads" --paired "${compression[@]}" \
    --report "$out/kraken2/$sample.report" \
    --output "$out/kraken2/$sample.kraken.output" "$r1" "$r2" \
    >"$out/logs/$sample.kraken.stdout.log" 2>"$out/logs/$sample.kraken.stderr.log"
  echo "[$sample] Bracken" >&2
  bracken -d "$kraken_db" -i "$out/kraken2/$sample.report" \
    -o "$out/bracken/$sample.bracken.tsv" -w "$out/bracken/$sample.bracken.report" \
    -r "$read_length" -l S >"$out/logs/$sample.bracken.log" 2>&1
  echo "[$sample] FASTQ preparation and HUMAnN" >&2
  combined="$out/humann_input/$sample.fastq"
  python3 "$script_dir/merge_paired_fastq.py" --r1 "$r1" --r2 "$r2" --output "$combined" \
    >"$out/logs/$sample.merge.log" 2>&1
  humann --input "$combined" --output "$out/humann/$sample" --threads "$threads" \
    --nucleotide-database "$nucleotide_db" --protein-database "$protein_db" \
    --metaphlan-options "--bowtie2db $metaphlan_db" \
    >"$out/logs/$sample.humann.stdout.log" 2>"$out/logs/$sample.humann.stderr.log"
  echo "[$sample] Done" >&2
done < <(tail -n +2 "$samples")
