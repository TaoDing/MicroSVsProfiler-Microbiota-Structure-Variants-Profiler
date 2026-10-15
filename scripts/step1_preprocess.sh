#!/usr/bin/env bash
# Step 1: KneadData QC and synchronized-seed seqtk downsampling.
# Usage: bash scripts/step1_preprocess.sh manifest.tsv HOST_DB OUTPUT_DIR [THREADS]
set -euo pipefail

if (( $# < 3 || $# > 4 )); then
  echo "Usage: $0 manifest.tsv HOST_BOWTIE2_DB OUTPUT_DIR [THREADS=20]" >&2
  exit 2
fi
manifest=$1
host_db=$2
out=$3
threads=${4:-20}
target_pairs=${TARGET_PAIRS:-10000000}
seed=${SEED:-100}

for program in kneaddata seqtk python3; do
  command -v "$program" >/dev/null || { echo "Missing executable: $program" >&2; exit 1; }
done
[[ -f "$manifest" ]] || { echo "Manifest not found: $manifest" >&2; exit 1; }
[[ -d "$host_db" ]] || { echo "Host reference directory not found: $host_db" >&2; exit 1; }
mkdir -p "$out/kneaddata" "$out/seqtk"

# TSV columns: sample_id,read1,read2 (with a header); paths may be absolute
# or relative to the directory containing the manifest.
manifest_dir=$(cd "$(dirname "$manifest")" && pwd)
while IFS=$'\t' read -r sample read1 read2 extra || [[ -n "${sample:-}" ]]; do
  [[ -z "${sample:-}" || "$sample" == 'sample_id' || "$sample" == \#* ]] && continue
  [[ -z "${extra:-}" ]] || { echo "Extra column for $sample" >&2; exit 1; }
  [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe sample ID: $sample" >&2; exit 1; }
  [[ -n "${read1:-}" && -n "${read2:-}" ]] || { echo "Missing read paths for $sample" >&2; exit 1; }
  [[ "$read1" = /* ]] || read1="$manifest_dir/$read1"
  [[ "$read2" = /* ]] || read2="$manifest_dir/$read2"
  [[ -f "$read1" && -f "$read2" ]] || { echo "Missing FASTQ for $sample" >&2; exit 1; }

  sample_dir="$out/kneaddata/$sample"
  mkdir -p "$sample_dir"
  echo "[$sample] KneadData QC..."
  kneaddata \
    -i "$read1" -i "$read2" \
    -db "$host_db" --reorder \
    -t "$threads" -p "$threads" \
    -o "$sample_dir" --log "$sample_dir/$sample.kneaddata.log"

  # KneadData 0.7.5 output naming depends on input FASTQ basename.
  # Require exactly one paired_1 and paired_2 file, and fail otherwise.
  mapfile -d '' mate1 < <(find "$sample_dir" -maxdepth 1 -type f \( -name '*_kneaddata_paired_1.fastq' -o -name '*_kneaddata_paired_1.fastq.gz' \) -print0)
  mapfile -d '' mate2 < <(find "$sample_dir" -maxdepth 1 -type f \( -name '*_kneaddata_paired_2.fastq' -o -name '*_kneaddata_paired_2.fastq.gz' \) -print0)
  (( ${#mate1[@]} == 1 && ${#mate2[@]} == 1 )) || {
    echo "[$sample] Expected exactly one paired FASTQ per mate; inspect $sample_dir" >&2
    exit 1
  }

  # Same seed and same number/order of paired records are essential.
  # Validate pairing with the accompanying Python validator before sampling.
  python3 "$(dirname "$0")/validate_fastq_pairs.py" "${mate1[0]}" "${mate2[0]}"
  echo "[$sample] seqtk: up to $target_pairs read pairs (seed $seed)..."
  seqtk sample -s "$seed" "${mate1[0]}" "$target_pairs" > "$out/seqtk/${sample}_R1.fastq"
  seqtk sample -s "$seed" "${mate2[0]}" "$target_pairs" > "$out/seqtk/${sample}_R2.fastq"
  python3 "$(dirname "$0")/validate_fastq_pairs.py" "$out/seqtk/${sample}_R1.fastq" "$out/seqtk/${sample}_R2.fastq"
  echo "[$sample] Done"
done < "$manifest"
