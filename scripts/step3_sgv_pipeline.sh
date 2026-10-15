#!/usr/bin/env bash
# Wrapper around the original SGV-Finder ICRA, PerFile and cohort-level CLI.
set -euo pipefail

usage() {
cat <<'HELP'
Usage:
  bash scripts/step3_sgv_pipeline.sh icra --sgv-src DIR --reads-prefix PREFIX --sample ID [--python PYTHON] [--out DIR]
  bash scripts/step3_sgv_pipeline.sh perfile --sgv-src DIR --jsdel FILE --sample ID [--read-length 150] [--python PYTHON] [--out DIR]
  bash scripts/step3_sgv_pipeline.sh call --sgv-src DIR --cohort NAME --input-glob 'PATTERN' [--min-samp-cutoff 75] [--python PYTHON] [--out DIR]

Defaults: --out results/step3, --python python, --read-length 150,
          --min-samp-cutoff 75

Important:
  * The legacy vendor Python programs must remain those supplied by SGV-Finder.
  * --reads-prefix is the input *prefix* expected by ICRA_cmd.py, not an R1 path.
  * --jsdel is the actual existing output path from ICRA; locate it in the
    ICRA output rather than assuming a particular filename.
  * --input-glob must be quoted and select the cohort's PerFile .jsdel files.
  * Test actual option compatibility against your installed SGV-Finder revision.
HELP
}
[[ $# -ge 1 ]] || { usage; exit 2; }
mode=$1; shift
src='' sample='' prefix='' jsdel='' cohort='' glob='' out='results/step3' py='python' readlen=150 cutoff=75
while (($#)); do
  case "$1" in
    --sgv-src|--sample|--reads-prefix|--jsdel|--cohort|--input-glob|--out|--python|--read-length|--min-samp-cutoff)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      k=$1; v=$2; shift 2
      case "$k" in
        --sgv-src) src=$v;; --sample) sample=$v;; --reads-prefix) prefix=$v;;
        --jsdel) jsdel=$v;; --cohort) cohort=$v;; --input-glob) glob=$v;;
        --out) out=$v;; --python) py=$v;; --read-length) readlen=$v;;
        --min-samp-cutoff) cutoff=$v;;
      esac;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2;;
  esac
done
[[ -n "$src" && -d "$src" ]] || { echo 'Missing or invalid --sgv-src' >&2; exit 2; }
[[ "$readlen" =~ ^[0-9]+$ && "$cutoff" =~ ^[0-9]+$ ]] || { echo 'Numerical flags must be integers' >&2; exit 2; }
mkdir -p "$out"
case "$mode" in
  icra)
    [[ -n "$sample" && -n "$prefix" ]] || { usage; exit 2; }
    [[ -f "$src/ICRA_cmd.py" ]] || { echo "Missing $src/ICRA_cmd.py" >&2; exit 2; }
    dest="$out/icra/$sample"; mkdir -p "$dest"
    "$py" "$src/ICRA_cmd.py" "$dest" "$prefix" --pe --use_theta --debug
    echo "ICRA completed. Inspect $dest for the actual .jsdel output before PerFile."
    ;;
  perfile)
    [[ -n "$sample" && -f "$jsdel" ]] || { echo 'Provide --sample and an existing --jsdel file' >&2; exit 2; }
    [[ -f "$src/SGVF_PerFile_cmd.py" ]] || { echo "Missing $src/SGVF_PerFile_cmd.py" >&2; exit 2; }
    dest="$out/perfile"; mkdir -p "$dest"
    "$py" "$src/SGVF_PerFile_cmd.py" "$jsdel" "$dest/${sample}.PerFile.jsdel" "$readlen" --x_coverage 0.01 --rate_param 10
    ;;
  call)
    [[ -n "$cohort" && -n "$glob" ]] || { usage; exit 2; }
    [[ "$cohort" =~ ^[a-zA-Z0-9_.-]+$ ]] || { echo 'Invalid cohort label' >&2; exit 2; }
    # The source document used SGVF_cmd_changed.py, which can be a local fork.
    [[ -f "$src/SGVF_cmd_changed.py" ]] || { echo "Missing $src/SGVF_cmd_changed.py; inspect your SGV-Finder revision" >&2; exit 2; }
    # Verify that the glob matches at least one file without evaluating arbitrary shell expressions.
    "$py" - "$glob" <<'PY'
import glob, sys
files = glob.glob(sys.argv[1])
if not files:
    sys.exit('No files matched --input-glob')
print('Matched %d cohort PerFile files' % len(files))
PY
    dest="$out/sv/$cohort"; mkdir -p "$dest"
    "$py" "$src/SGVF_cmd_changed.py" "$glob" "$dest/${cohort}.dsgv.csv" "$dest/${cohort}.vsgv.csv" \
      --min_samp_cutoff "$cutoff" --x_coverage 0.01 --rate_param 10 \
      --browser_path "$dest/browser" --csv_output --byorig
    ;;
  *) echo "Unknown subcommand: $mode" >&2; usage; exit 2;;
esac
