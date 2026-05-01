#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[$(date -Is)] $*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

usage(){
  cat <<'USAGE'
ONT amplicon hybrid reference benchmark, per barcode folder

Fast strategy:
  1) Trim and filter reads once per sample.
  2) Screen the full multi-reference panel once with minimap2.
  3) Shortlist top references from the screen.
  4) Benchmark shortlisted references independently.
  5) Run consensus only on shortlisted candidates.

Required:
  --refs refs.fasta
  --primers-xlsx BarcodedPrimers.xlsx
  --barcode-map barcodes.csv
  --fastq-pass /path/to/fastq_pass
  --out outdir

Optional:
  --preset NAME                  Preset profile. Supported: piranha-vp1
  --piranha-vp1-mode             Alias for --preset piranha-vp1
  --threads N                     Default all available CPU threads
  --min-len 800                   Default 800
  --min-read-length N             Alias for --min-len
  --max-len 1200                  Default 1200
  --max-read-length N             Alias for --max-len
  --min-q 10                      Default 10
  --min-cov 50                    Default 50
  --include-unclassified 0|1      Default 0
  --medaka-model auto|MODEL       Default auto
  --medaka-resolve-reads N        Default 500
  --primer-trim-mode MODE         exact|fixed|none. Default exact
  --primer-length N               Fixed-end trim length for MODE=fixed. Default 30
  --min-map-quality N             Minimum MAPQ used for counts/consensus masking. Default 0
  --min-base-quality N            Minimum base quality used for depth masking. Default 0
  --minimap2-options "k=v ..."    minimap2 overrides used for screening and benchmarking. Default x=map-ont
  --min-aln-block N               Minimum alignment block length to count a hit. Default 0
  --exclude-supplementary 0|1     Default 1
  --screen-topN N                 Default 20
  --screen-min-mapped N           Default 20
  --rank-map-topN N               Default 10. Use 0 to run consensus on all shortlisted refs
  --consensus-all-passing 0|1     Default 0
  --attempt-consensus-any-mapped 0|1  Default 0. Attempt consensus for every shortlisted ref with >0 benchmark mapped reads, ignoring breadth/read thresholds.
  --min-breadth-for-consensus F   Default 0.40
  --min-mapped-for-consensus N    Default 50
  --min-read-depth N              Alias for --min-mapped-for-consensus
  --min-read-pcent F              Minimum percent of sample required for consensus attempt. Default 0
  --min-cons-len N                Default 800
  --max-n-frac F                  Default 0.30
  --final-topK N                  Default 3
  --keep-all-bams 0|1             Default 0
  --run-report-html PATH           Optional MinKNOW run report html to include in HTML report
  --run-name NAME                  Optional report title/name

Outputs:
  outdir/screen_per_reference.tsv
  outdir/benchmark_per_reference.tsv
  outdir/summary.tsv
  outdir/overall_reference_leaderboard.tsv
  outdir/all_best_consensus.fasta
  outdir/top_hits_report.tsv
  outdir/report.html
  outdir/results_bundle.zip
USAGE
}

need_cmd(){ command -v "$1" >/dev/null 2>&1; }

detect_threads(){
  local n=""
  if need_cmd nproc; then n="$(nproc 2>/dev/null || true)"; fi
  if [[ -z "$n" ]] && need_cmd getconf; then n="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"; fi
  if [[ -z "$n" ]] && need_cmd python3; then
    n="$(python3 - <<'PY'
import os
print(os.cpu_count() or 1)
PY
)"
  fi
  if [[ ! "$n" =~ ^[0-9]+$ || "$n" -lt 1 ]]; then n=1; fi
  echo "$n"
}

apply_piranha_vp1_preset(){
  MIN_LEN=1000
  MAX_LEN=1300
  MIN_Q=0
  MEDAKA_MODEL="r1041_e82_400bps_hac_variant_v4.3.0"
  TRIM_MODE="fixed"
  PRIMER_LENGTH=30
  MIN_MAPQ=0
  MIN_BASEQ=9
  MINIMAP2_OPTIONS="x=map-ont"
  MIN_ALN_BLOCK=600
  MIN_BREADTH_FOR_CONS=0
  MIN_MAPPED_FOR_CONS=50
  MIN_MAPPED_PCENT_FOR_CONS=0
  MIN_CONS_LEN=0
  MAX_N_FRAC=1
}

build_minimap2_args(){
  local out_name="$1"
  local option_string="${2:-}"
  local default_preset="${3:-map-ont}"
  local token flag value
  local have_x=0
  local -n out_arr="$out_name"
  out_arr=(-a)
  if [[ -n "$option_string" ]]; then
    read -r -a _mm2_tokens <<< "$option_string"
    for token in "${_mm2_tokens[@]}"; do
      token="${token#,}"
      token="${token%,}"
      [[ -n "$token" ]] || continue
      [[ "$token" == *=* ]] || die "minimap2 option items must be flag=value, got: $token"
      flag="${token%%=*}"
      value="${token#*=}"
      [[ "$flag" =~ ^[A-Za-z]$ ]] || die "Unsupported minimap2 option key: $flag"
      out_arr+=("-$flag" "$value")
      [[ "$flag" == "x" ]] && have_x=1
    done
  fi
  if [[ "$have_x" -eq 0 ]]; then
    out_arr+=(-x "$default_preset")
  fi
  out_arr+=(--secondary=no -N 1)
}

filter_bam_by_min_aln_block(){
  local in_bam="$1"
  local out_bam="$2"
  local min_aln="$3"
  samtools view -h "$in_bam"     | awk -v mina="$min_aln" '
        function alnlen(c,    total,seg,n,op) {
          total=0
          while (match(c, /[0-9]+[MIDNSHP=X]/)) {
            seg=substr(c, RSTART, RLENGTH)
            n=seg+0
            op=substr(seg, length(seg), 1)
            if (op ~ /[MIDN=X]/) total += n
            c=substr(c, RSTART + RLENGTH)
          }
          return total
        }
        /^@/ { print; next }
        { if (alnlen($6) > mina) print }
      '     | samtools view -b -o "$out_bam" -
}

THREADS="$(detect_threads)"
MIN_LEN=800
MAX_LEN=1200
MIN_Q=10
MIN_COV=50
INCLUDE_UNCLASSIFIED=0
MEDAKA_MODEL="auto"
MEDAKA_RESOLVE_READS=500
TRIM_MODE="exact"
PRIMER_LENGTH=30
MIN_MAPQ=0
MIN_BASEQ=0
MINIMAP2_OPTIONS=""
MIN_ALN_BLOCK=0
EXCLUDE_SUPP=1
SCREEN_TOPN=20
SCREEN_MIN_MAPPED=20
RANK_MAP_TOPN=10
CONSENSUS_ALL_PASSING=0
ATTEMPT_CONSENSUS_ANY_MAPPED=0
MIN_BREADTH_FOR_CONS=0.40
MIN_MAPPED_FOR_CONS=50
MIN_MAPPED_PCENT_FOR_CONS=0
MIN_CONS_LEN=800
MAX_N_FRAC=0.30
FINAL_TOPK=3
KEEP_ALL_BAMS=0
RUN_REPORT_HTML=""
RUN_NAME=""

PRESET=""
REFS=""
PRIMERS_XLSX=""
BARCODE_MAP=""
FASTQ_PASS=""
OUT=""

ARGS=("$@")
for ((i=0; i<${#ARGS[@]}; i++)); do
  case "${ARGS[$i]}" in
    --preset)
      (( i + 1 < ${#ARGS[@]} )) || die "--preset requires a value"
      PRESET="${ARGS[$((i+1))]}"
      ;;
    --piranha-vp1-mode)
      PRESET="piranha-vp1"
      ;;
  esac
done

case "$PRESET" in
  "") ;;
  piranha-vp1) apply_piranha_vp1_preset ;;
  *) die "Unsupported --preset value: $PRESET" ;;
esac

while [[ $# -gt 0 ]]; do
  case "$1" in
    --preset) PRESET="$2"; shift 2;;
    --piranha-vp1-mode) PRESET="piranha-vp1"; shift;;
    --refs) REFS="$2"; shift 2;;
    --primers-xlsx) PRIMERS_XLSX="$2"; shift 2;;
    --barcode-map) BARCODE_MAP="$2"; shift 2;;
    --fastq-pass) FASTQ_PASS="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --min-len) MIN_LEN="$2"; shift 2;;
    --min-read-length) MIN_LEN="$2"; shift 2;;
    --max-len) MAX_LEN="$2"; shift 2;;
    --max-read-length) MAX_LEN="$2"; shift 2;;
    --min-q) MIN_Q="$2"; shift 2;;
    --min-cov) MIN_COV="$2"; shift 2;;
    --include-unclassified) INCLUDE_UNCLASSIFIED="$2"; shift 2;;
    --medaka-model) MEDAKA_MODEL="$2"; shift 2;;
    --medaka-resolve-reads) MEDAKA_RESOLVE_READS="$2"; shift 2;;
    --primer-trim-mode) TRIM_MODE="$2"; shift 2;;
    --primer-length) PRIMER_LENGTH="$2"; shift 2;;
    --min-map-quality) MIN_MAPQ="$2"; shift 2;;
    --min-base-quality) MIN_BASEQ="$2"; shift 2;;
    --minimap2-options) MINIMAP2_OPTIONS="$2"; shift 2;;
    --min-aln-block) MIN_ALN_BLOCK="$2"; shift 2;;
    --exclude-supplementary) EXCLUDE_SUPP="$2"; shift 2;;
    --screen-topN) SCREEN_TOPN="$2"; shift 2;;
    --screen-min-mapped) SCREEN_MIN_MAPPED="$2"; shift 2;;
    --rank-map-topN) RANK_MAP_TOPN="$2"; shift 2;;
    --consensus-all-passing) CONSENSUS_ALL_PASSING="$2"; shift 2;;
    --attempt-consensus-any-mapped) ATTEMPT_CONSENSUS_ANY_MAPPED="$2"; shift 2;;
    --min-breadth-for-consensus) MIN_BREADTH_FOR_CONS="$2"; shift 2;;
    --min-mapped-for-consensus) MIN_MAPPED_FOR_CONS="$2"; shift 2;;
    --min-read-depth) MIN_MAPPED_FOR_CONS="$2"; shift 2;;
    --min-read-pcent) MIN_MAPPED_PCENT_FOR_CONS="$2"; shift 2;;
    --min-cons-len) MIN_CONS_LEN="$2"; shift 2;;
    --max-n-frac) MAX_N_FRAC="$2"; shift 2;;
    --final-topK) FINAL_TOPK="$2"; shift 2;;
    --keep-all-bams) KEEP_ALL_BAMS="$2"; shift 2;;
    --run-report-html) RUN_REPORT_HTML="$2"; shift 2;;
    --run-name) RUN_NAME="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$REFS" && -f "$REFS" ]] || { usage; die "--refs is required and must exist"; }
[[ -n "$PRIMERS_XLSX" && -f "$PRIMERS_XLSX" ]] || { usage; die "--primers-xlsx is required and must exist"; }
[[ -n "$BARCODE_MAP" && -f "$BARCODE_MAP" ]] || { usage; die "--barcode-map is required and must exist"; }
[[ -n "$FASTQ_PASS" && -d "$FASTQ_PASS" ]] || { usage; die "--fastq-pass is required and must exist"; }
[[ -n "$OUT" ]] || { usage; die "--out is required"; }
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 ]] || die "--threads must be an integer >= 1"
[[ "$TRIM_MODE" =~ ^(exact|fixed|none)$ ]] || die "--primer-trim-mode must be one of: exact, fixed, none"
[[ "$PRIMER_LENGTH" =~ ^[0-9]+$ ]] || die "--primer-length must be an integer >= 0"
[[ "$MIN_MAPQ" =~ ^[0-9]+$ ]] || die "--min-map-quality must be an integer >= 0"
[[ "$MIN_BASEQ" =~ ^[0-9]+$ ]] || die "--min-base-quality must be an integer >= 0"
[[ "$MIN_ALN_BLOCK" =~ ^[0-9]+$ ]] || die "--min-aln-block must be an integer >= 0"
[[ "$MIN_MAPPED_PCENT_FOR_CONS" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--min-read-pcent must be a number >= 0"

build_minimap2_args MINIMAP2_ARGS "$MINIMAP2_OPTIONS" "map-ont"

req=(python3 minimap2 samtools cutadapt filtlong medaka medaka_consensus)
missing=0
for c in "${req[@]}"; do
  if ! need_cmd "$c"; then
    log "Missing dependency: $c"
    missing=1
  fi
done
if [[ "$missing" -eq 1 ]]; then
  exit 2
fi

mkdir -p "$OUT"/{00_logs,01_merge,02_trim,03_filter,04_screen,05_benchmark,06_consensus,per_sample,primers,ref,ref_cache}

{
  echo "date=$(date -Is)"
  echo "preset=$PRESET"
  echo "threads=$THREADS"
  echo "min_len=$MIN_LEN"
  echo "max_len=$MAX_LEN"
  echo "min_q=$MIN_Q"
  echo "min_cov=$MIN_COV"
  echo "screen_topN=$SCREEN_TOPN"
  echo "screen_min_mapped=$SCREEN_MIN_MAPPED"
  echo "rank_map_topN=$RANK_MAP_TOPN"
  echo "consensus_all_passing=$CONSENSUS_ALL_PASSING"
  echo "attempt_consensus_any_mapped=$ATTEMPT_CONSENSUS_ANY_MAPPED"
  echo "min_breadth_for_consensus=$MIN_BREADTH_FOR_CONS"
  echo "min_mapped_for_consensus=$MIN_MAPPED_FOR_CONS"
  echo "min_read_pcent_for_consensus=$MIN_MAPPED_PCENT_FOR_CONS"
  echo "min_cons_len=$MIN_CONS_LEN"
  echo "max_n_frac=$MAX_N_FRAC"
  echo "final_topK=$FINAL_TOPK"
  echo "primer_trim_mode=$TRIM_MODE"
  echo "primer_length=$PRIMER_LENGTH"
  echo "min_map_quality=$MIN_MAPQ"
  echo "min_base_quality=$MIN_BASEQ"
  echo "min_aln_block=$MIN_ALN_BLOCK"
  echo "minimap2_options=${MINIMAP2_OPTIONS:-x=map-ont}"
  echo "exclude_supplementary=$EXCLUDE_SUPP"
  echo "keep_all_bams=$KEEP_ALL_BAMS"
  echo "python=$(python3 --version 2>/dev/null || true)"
  echo "minimap2=$(minimap2 --version 2>/dev/null || true)"
  echo "samtools=$(samtools --version 2>/dev/null | head -n1 || true)"
  echo "cutadapt=$(cutadapt --version 2>/dev/null || true)"
  echo "filtlong=$(filtlong --version 2>/dev/null | head -n1 || true)"
  echo "medaka=$(medaka --version 2>/dev/null | head -n1 || true)"
} > "$OUT/00_logs/versions.txt"

BARCODE_TSV="$OUT/00_logs/barcode_map.tsv"
python3 - "$BARCODE_MAP" "$BARCODE_TSV" <<'PY'
import sys, csv, re
inp, outp = sys.argv[1], sys.argv[2]
rows=[]
seen=set()
def normalize_barcode(v):
    s=str(v).strip()
    if not s:
        return ""
    sl=s.lower()
    if sl in {"unclassified", "uncl"}:
        return "unclassified"
    m=re.search(r"(\d+)", s)
    if m:
        return f"barcode{int(m.group(1)):02d}"
    return s
with open(inp, newline="") as f:
    r=csv.DictReader(f)
    if not r.fieldnames or "barcode" not in r.fieldnames or "sample" not in r.fieldnames:
        raise SystemExit("barcodes.csv must have columns: barcode,sample")
    for row in r:
        b=normalize_barcode(row["barcode"])
        s=str(row["sample"]).strip()
        if b and b not in seen:
            rows.append((b, s))
            seen.add(b)
with open(outp, "w", newline="") as g:
    g.write("barcode\tsample\n")
    for b,s in rows:
        g.write(f"{b}\t{s}\n")
PY

PRIMERS_BC="$OUT/primers/primers_by_bc.tsv"
python3 - "$PRIMERS_XLSX" "$PRIMERS_BC" <<'PY'
import sys, pandas as pd
xlsx, out_tsv = sys.argv[1], sys.argv[2]
df = pd.read_excel(xlsx, sheet_name=0)
required = {"Name","Orientation","Sequence"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Missing columns in Excel: {sorted(missing)}")
df = df.dropna(subset=["Name","Orientation","Sequence"])
df["Orientation"] = df["Orientation"].astype(str).str.strip().str.lower()
df["bc"] = df["Name"].astype(str).str.extract(r"__(BC\d+)$")[0]
df = df.dropna(subset=["bc"])
df["bc_num"] = df["bc"].str.replace("BC","", regex=False).astype(int)
df["bc"] = df["bc_num"].map(lambda x: f"BC{x:02d}")
out=[]
for bc, g in df.groupby("bc"):
    f = g.loc[g["Orientation"]=="forward","Sequence"].dropna().unique().tolist()
    r = g.loc[g["Orientation"]=="reverse","Sequence"].dropna().unique().tolist()
    if len(f)!=1 or len(r)!=1:
        raise SystemExit(f"{bc} needs exactly 1 forward and 1 reverse primer. forward={len(f)} reverse={len(r)}")
    out.append((bc, f[0].strip(), r[0].strip()))
out = sorted(out, key=lambda x: int(x[0].replace("BC","")))
with open(out_tsv,"w") as w:
    w.write("bc\tfwd_seq\trev_seq\n")
    for bc,f,r in out:
        w.write(f"{bc}\t{f}\t{r}\n")
PY

revcomp(){ echo "$1" | tr 'ACGTacgtNn' 'TGCAtgcaNn' | rev; }

log "Indexing full multi-FASTA reference once"
samtools faidx "$REFS"
minimap2 -d "$OUT/ref_cache/refs.full.mmi" "$REFS" >/dev/null 2>&1 || true

awk '
  BEGIN{OFS="\t"}
  /^>/{h=$0; sub(/^>/,"",h); split(h,a," "); id=a[1]; print id, h}
' "$REFS" > "$OUT/ref/ref_headers.tsv"

REF_META="$OUT/ref_cache/ref_meta.tsv"
python3 - "$REFS" "$REF_META" <<'PY'
import sys, re
inp, outp = sys.argv[1], sys.argv[2]
def parse_header_meta(header):
    ddns = ""
    cluster = ""
    m = re.search(r"(?:^|\s)ddns_group=([^\s]+)", header)
    if m: ddns = m.group(1)
    m = re.search(r"(?:^|\s)cluster=([^\s]+)", header)
    if m: cluster = m.group(1)
    return ddns, cluster
name=None; header=None; buf=[]
rows=[]
with open(inp) as f:
    for line in f:
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                seq="".join(buf).upper()
                ddns, cluster = parse_header_meta(header)
                rows.append((name, len(seq), header, ddns or "Other", cluster or name))
            header=line[1:]
            name=header.split()[0]
            buf=[]
        else:
            buf.append(line.strip())
    if name is not None:
        seq="".join(buf).upper()
        ddns, cluster = parse_header_meta(header)
        rows.append((name, len(seq), header, ddns or "Other", cluster or name))
with open(outp, "w") as w:
    w.write("ref_id\tref_len\theader\tddns_group\tcluster\n")
    for r in rows:
        w.write("\t".join(map(str, r)) + "\n")
PY

MEDAKA_MODELS_LIST="$(medaka tools list_models 2>/dev/null | sed 's/^Available:[[:space:]]*//' | tr ',' '\n' | tr -d ' ' || true)"
RUN_MEDAKA_MODEL=""
append_model_to_versions(){
  local m="$1"
  if [[ -f "$OUT/00_logs/versions.txt" ]] && ! grep -q '^medaka_model=' "$OUT/00_logs/versions.txt" 2>/dev/null; then
    echo "medaka_model=$m" >> "$OUT/00_logs/versions.txt"
  fi
}
pick_medaka_fallback(){
  local m=""
  m="$(echo "$MEDAKA_MODELS_LIST" | grep -E '^r1041_.*_400bps_.*hac' | head -n1 || true)"
  if [[ -z "$m" ]]; then m="$(echo "$MEDAKA_MODELS_LIST" | grep -E '^r1041_.*_400bps_.*sup' | head -n1 || true)"; fi
  if [[ -z "$m" ]]; then m="$(echo "$MEDAKA_MODELS_LIST" | grep -E '^r941_.*_hac' | head -n1 || true)"; fi
  if [[ -z "$m" ]]; then m="$(echo "$MEDAKA_MODELS_LIST" | head -n1 || true)"; fi
  echo "$m"
}
make_model_sample_fastq(){
  local in_gz="$1"
  local out_fastq="$2"
  local nreads="$3"
  python3 - "$in_gz" "$out_fastq" "$nreads" <<'PY'
import sys, gzip
inp, outp, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
def open_in(p):
    if p.endswith('.gz'): return gzip.open(p, 'rt', errors='replace')
    return open(p, 'rt', errors='replace')
reads = 0
with open_in(inp) as f, open(outp, 'w') as w:
    while reads < n:
        h = f.readline()
        if not h: break
        s = f.readline(); p = f.readline(); q = f.readline()
        if not q: break
        w.write(h); w.write(s); w.write(p); w.write(q)
        reads += 1
print(reads)
PY
}
resolve_medaka_model_once(){
  local fastq_gz="$1"
  [[ -s "$fastq_gz" ]] || return 1
  local tmp="$OUT/00_logs/medaka_model_sample.fastq"
  local n
  n="$(make_model_sample_fastq "$fastq_gz" "$tmp" "$MEDAKA_RESOLVE_READS" 2>/dev/null | tail -n1 || true)"
  [[ -n "$n" && "$n" -gt 0 ]] || return 1
  local resolved
  resolved="$(medaka tools resolve_model --auto_model consensus "$tmp" 2>>"$OUT/00_logs/medaka_model_resolve.log" | tail -n1 | tr -d '[:space:]' || true)"
  [[ -n "$resolved" ]] || return 1
  if [[ -f "$resolved" ]]; then
    local b
    b="$(basename "$resolved")"
    b="${b%_model_pt.tar.gz}"
    b="${b%_model.tar.gz}"
    resolved="$b"
  fi
  RUN_MEDAKA_MODEL="$resolved"
  echo "$RUN_MEDAKA_MODEL" > "$OUT/00_logs/medaka_model.txt"
  append_model_to_versions "$RUN_MEDAKA_MODEL"
  return 0
}

if [[ -n "$MEDAKA_MODEL" && "$MEDAKA_MODEL" != "auto" ]]; then
  RUN_MEDAKA_MODEL="$MEDAKA_MODEL"
  echo "$RUN_MEDAKA_MODEL" > "$OUT/00_logs/medaka_model.txt"
  append_model_to_versions "$RUN_MEDAKA_MODEL"
fi

SCREEN_MASTER="$OUT/screen_per_reference.tsv"
BENCH_MASTER="$OUT/benchmark_per_reference.tsv"
CONS_MASTER="$OUT/00_logs/per_sample_consensus_metrics.tsv"
echo -e "sample_id\tsample_name\tbarcode_folder\tbc\treads_raw\treads_trimmed\treads_filtered\tref_id\tddns_group\tcluster\tref_len\tscreen_mapped_reads\tscreen_pct_mapped" > "$SCREEN_MASTER"
echo -e "sample_id\tsample_name\tbarcode_folder\tbc\treads_raw\treads_trimmed\treads_filtered\tref_id\tddns_group\tcluster\tref_len\tmapped_reads\tpct_mapped\tbreadth\tmean_depth\tmin_depth\tmap_bam\tdepth_tsv" > "$BENCH_MASTER"
echo -e "sample_id\tref_id\tconsensus_attempted\tconsensus_status_raw\tconsensus_fasta\tconsensus_len\tn_count\tcallable_bases\tmutations\tmut_rate\tmedaka_bam" > "$CONS_MASTER"

FILT_FLAG=260
if [[ "$EXCLUDE_SUPP" -eq 1 ]]; then FILT_FLAG=2308; fi
VIEW_FILTER_ARGS=(-F "$FILT_FLAG")
if [[ "$MIN_MAPQ" -gt 0 ]]; then VIEW_FILTER_ARGS+=(-q "$MIN_MAPQ"); fi
DEPTH_ARGS=(-aa)
if [[ "$MIN_MAPQ" -gt 0 ]]; then DEPTH_ARGS+=(-q "$MIN_MAPQ"); fi
if [[ "$MIN_BASEQ" -gt 0 ]]; then DEPTH_ARGS+=(-Q "$MIN_BASEQ"); fi

mapfile -t folders < <(awk -F'\t' 'NR>1 && $1!="" && !seen[$1]++ {print $1}' "$BARCODE_TSV")
if [[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]]; then
  if ! printf '%s\n' "${folders[@]}" | grep -Fxq "unclassified"; then folders+=("unclassified"); fi
fi
[[ ${#folders[@]} -gt 0 ]] || die "No barcode entries found"

for barcode_folder in "${folders[@]}"; do
  [[ -n "$barcode_folder" ]] || continue
  d="$FASTQ_PASS/$barcode_folder"
  [[ -d "$d" ]] || { log "missing barcode folder, skipping: $barcode_folder"; continue; }

  bc="NA"
  if [[ "$barcode_folder" =~ ^barcode([0-9]+)$ ]]; then
    bc=$(printf "BC%02d" "$((10#${BASH_REMATCH[1]}))")
  elif [[ "$barcode_folder" == "unclassified" ]]; then
    bc="UNCL"
  fi

  sample_name="$(awk -F'\t' -v b="$barcode_folder" 'NR>1 && $1==b {print $2; exit}' "$BARCODE_TSV" || true)"
  [[ -n "$sample_name" ]] || sample_name="$barcode_folder"
  safe_sample="$(echo "$sample_name" | sed 's/[^A-Za-z0-9._-]/_/g')"
  sample_id="${safe_sample}__${barcode_folder}"
  log "Processing $sample_id"

  shopt -s nullglob
  fqs=( "$d"/*.fastq.gz "$d"/*.fq.gz "$d"/*.fastq "$d"/*.fq )
  [[ ${#fqs[@]} -gt 0 ]] || { log "no FASTQ files, skipping"; continue; }

  merged="$OUT/01_merge/${sample_id}.fastq.gz"
  if [[ ! -s "$merged" ]]; then
    { for f in "${fqs[@]}"; do if [[ "$f" == *.gz ]]; then zcat "$f"; else cat "$f"; fi; done; } | gzip -c > "$merged"
  fi
  reads_raw=$(zcat "$merged" | awk 'END{print NR/4}')

  trimmed="$merged"
  reads_trimmed="$reads_raw"

  if [[ "$TRIM_MODE" == "exact" && "$bc" != "UNCL" ]]; then
    row="$(awk -F'	' -v bc="$bc" 'NR>1 && $1==bc {print; exit}' "$PRIMERS_BC" || true)"
    [[ -n "$row" ]] || { log "no primers for $bc, skipping"; continue; }
    fwd="$(echo "$row" | cut -f2)"
    rev="$(echo "$row" | cut -f3)"
    fwd_rc="$(revcomp "$fwd")"
    rev_rc="$(revcomp "$rev")"
    t1="$OUT/02_trim/${sample_id}.pass1.fastq.gz"
    u1="$OUT/02_trim/${sample_id}.pass1.untrim.fastq.gz"
    cutadapt -j "$THREADS" -e 0.15 --overlap 20 -g "$fwd" -a "$rev_rc" --untrimmed-output "$u1" -o "$t1" "$merged" > "$OUT/00_logs/${sample_id}.cutadapt.pass1.log" 2>&1
    t2="$OUT/02_trim/${sample_id}.pass2.fastq.gz"
    cutadapt -j "$THREADS" -e 0.15 --overlap 20 -g "$rev" -a "$fwd_rc" --discard-untrimmed -o "$t2" "$u1" > "$OUT/00_logs/${sample_id}.cutadapt.pass2.log" 2>&1
    trimmed="$OUT/02_trim/${sample_id}.trim.fastq.gz"
    { zcat "$t1" 2>/dev/null || true; zcat "$t2" 2>/dev/null || true; } | gzip -c > "$trimmed"
    reads_trimmed=$(zcat "$trimmed" | awk 'END{print NR/4}')
    if [[ "$reads_trimmed" -eq 0 ]]; then trimmed="$merged"; reads_trimmed="$reads_raw"; fi
  fi

  filter_input="$trimmed"
  if [[ "$TRIM_MODE" != "exact" ]]; then
    filter_input="$merged"
    reads_trimmed="$reads_raw"
  fi

  filtered="$OUT/03_filter/${sample_id}.filt.fastq.gz"
  if awk "BEGIN{exit !("$MIN_Q" > 0)}"; then
    filtlong --min_length "$MIN_LEN" --max_length "$MAX_LEN" --min_mean_q "$MIN_Q" "$filter_input" | gzip -c > "$filtered"
  else
    filtlong --min_length "$MIN_LEN" --max_length "$MAX_LEN" "$filter_input" | gzip -c > "$filtered"
  fi
  reads_filtered=$(zcat "$filtered" | awk 'END{print NR/4}')
  [[ "$reads_filtered" -gt 0 ]] || { log "No reads after filtering for $sample_id"; continue; }

  screen_reads="$filtered"
  benchmark_reads="$filtered"
  if [[ "$TRIM_MODE" == "fixed" && "$PRIMER_LENGTH" -gt 0 ]]; then
    fixed_trimmed="$OUT/03_filter/${sample_id}.filt.fixedtrim.fastq.gz"
    cutadapt -j "$THREADS" -u "$PRIMER_LENGTH" -u "-$PRIMER_LENGTH" -o "$fixed_trimmed" "$filtered" > "$OUT/00_logs/${sample_id}.cutadapt.fixedtrim.log" 2>&1 || true
    fixed_reads=$(zcat "$fixed_trimmed" 2>/dev/null | awk 'END{print NR/4}')
    if [[ "$fixed_reads" -gt 0 ]]; then
      benchmark_reads="$fixed_trimmed"
      reads_trimmed="$fixed_reads"
    else
      reads_trimmed="$reads_filtered"
    fi
  fi

  if [[ -z "$RUN_MEDAKA_MODEL" && ( -z "$MEDAKA_MODEL" || "$MEDAKA_MODEL" == "auto" ) ]]; then
    if ! resolve_medaka_model_once "$merged"; then
      RUN_MEDAKA_MODEL="$(pick_medaka_fallback)"
      [[ -n "$RUN_MEDAKA_MODEL" ]] || die "Could not select a medaka model"
      echo "$RUN_MEDAKA_MODEL" > "$OUT/00_logs/medaka_model.txt"
      append_model_to_versions "$RUN_MEDAKA_MODEL"
    fi
  fi

  sample_dir="$OUT/per_sample/$sample_id"
  mkdir -p "$sample_dir"/{screen,benchmark,consensus}

  screen_bam="$sample_dir/screen/allrefs.bam"
  screen_idx="$sample_dir/screen/idxstats.tsv"
  raw_screen_bam="$sample_dir/screen/allrefs.raw.bam"
  minimap2 -t "$THREADS" "${MINIMAP2_ARGS[@]}" "$OUT/ref_cache/refs.full.mmi" <(zcat "$screen_reads")     | samtools view -b "${VIEW_FILTER_ARGS[@]}"     | samtools sort -@ "$THREADS" -o "$raw_screen_bam"
  if [[ "$MIN_ALN_BLOCK" -gt 0 ]]; then
    filter_bam_by_min_aln_block "$raw_screen_bam" "$screen_bam" "$MIN_ALN_BLOCK"
    rm -f "$raw_screen_bam"
  else
    mv "$raw_screen_bam" "$screen_bam"
  fi
  samtools index "$screen_bam"
  samtools idxstats "$screen_bam" > "$screen_idx"
  python3 - "$screen_idx" "$sample_dir/idxstats.full.tsv" "$reads_filtered" <<'PY'
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
reads_filtered = int(sys.argv[3])

rows = []
mapped_total = 0

with src.open() as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        ref_id = parts[0]
        if ref_id == "*":
            continue
        mapped = int(parts[2])
        mapped_total += mapped
        rows.append(parts[:4])

reconstructed_unmapped = max(0, reads_filtered - mapped_total)

with dst.open("w") as out:
    for parts in rows:
        out.write("\t".join(parts) + "\n")
    out.write(f"*\t0\t0\t{reconstructed_unmapped}\n")
PY

  shortlist="$sample_dir/screen/shortlist.tsv"
  python3 - "$screen_idx" "$REF_META" "$reads_filtered" "$SCREEN_TOPN" "$SCREEN_MIN_MAPPED" "$SCREEN_MASTER" "$sample_id" "$sample_name" "$barcode_folder" "$bc" "$reads_raw" "$reads_trimmed" <<'PY' > "$shortlist"
import sys, csv
idx, meta, reads_filtered, topn, minm, master, sample_id, sample_name, barcode_folder, bc, reads_raw, reads_trimmed = sys.argv[1:]
reads_filtered = int(reads_filtered); topn = int(topn); minm = int(minm)
meta_d = {}
with open(meta) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        meta_d[row['ref_id']] = row
rows = []
with open(idx) as f:
    for line in f:
        ref_id, ref_len, mapped, unmapped = line.rstrip('\n').split('\t')
        if ref_id == '*' or ref_id not in meta_d:
            continue
        mapped = int(mapped); ref_len = int(ref_len)
        pct = 0.0 if reads_filtered == 0 else 100.0 * mapped / reads_filtered
        m = meta_d[ref_id]
        rows.append({
            'ref_id': ref_id,
            'ddns_group': m['ddns_group'],
            'cluster': m['cluster'],
            'ref_len': ref_len,
            'screen_mapped_reads': mapped,
            'screen_pct_mapped': pct,
        })
rows.sort(key=lambda x: (x['screen_mapped_reads'], x['screen_pct_mapped']), reverse=True)
chosen = [r for r in rows if r['screen_mapped_reads'] >= minm]
if topn > 0:
    chosen = chosen[:topn]
if not chosen and rows and rows[0]['screen_mapped_reads'] > 0:
    chosen = [rows[0]]
with open(master, 'a') as w:
    for r in rows:
        w.write(f"{sample_id}\t{sample_name}\t{barcode_folder}\t{bc}\t{reads_raw}\t{reads_trimmed}\t{reads_filtered}\t{r['ref_id']}\t{r['ddns_group']}\t{r['cluster']}\t{r['ref_len']}\t{r['screen_mapped_reads']}\t{r['screen_pct_mapped']:.2f}\n")
for r in chosen:
    sys.stdout.write(r['ref_id'] + '\n')
PY

  sample_map_tsv="$sample_dir/benchmark/map_metrics.tsv"
  sample_cons_tsv="$sample_dir/benchmark/consensus_metrics.tsv"
  echo -e "sample_id\tsample_name\tbarcode_folder\tbc\treads_raw\treads_trimmed\treads_filtered\tref_id\tddns_group\tcluster\tref_len\tmapped_reads\tpct_mapped\tbreadth\tmean_depth\tmin_depth\tmap_bam\tdepth_tsv" > "$sample_map_tsv"
  echo -e "sample_id\tref_id\tconsensus_attempted\tconsensus_status_raw\tconsensus_fasta\tconsensus_len\tn_count\tcallable_bases\tmutations\tmut_rate\tmedaka_bam" > "$sample_cons_tsv"

  while IFS= read -r ref_id; do
    [[ -n "$ref_id" ]] || continue
    meta_row="$(awk -F'\t' -v r="$ref_id" 'NR>1 && $1==r {print; exit}' "$REF_META" || true)"
    [[ -n "$meta_row" ]] || continue
    ref_len="$(echo "$meta_row" | cut -f2)"
    ddns_group="$(echo "$meta_row" | cut -f4)"
    cluster="$(echo "$meta_row" | cut -f5)"
    safe_name="$(echo "$ref_id" | sed 's/[^A-Za-z0-9._-]/_/g')"
    ref_fa="$OUT/ref_cache/${safe_name}.fa"
    if [[ ! -s "$ref_fa" ]]; then
      samtools faidx "$REFS" "$ref_id" > "$ref_fa"
      samtools faidx "$ref_fa"
      minimap2 -d "${ref_fa}.mmi" "$ref_fa" >/dev/null 2>&1 || true
    fi

    bam="$sample_dir/benchmark/${safe_name}.bam"
    depthfile="$sample_dir/benchmark/${safe_name}.depth.tsv"
    raw_bam="$sample_dir/benchmark/${safe_name}.raw.bam"
    minimap2 -t "$THREADS" "${MINIMAP2_ARGS[@]}" "${ref_fa}.mmi" <(zcat "$benchmark_reads")       | samtools view -b "${VIEW_FILTER_ARGS[@]}"       | samtools sort -@ "$THREADS" -o "$raw_bam"
    if [[ "$MIN_ALN_BLOCK" -gt 0 ]]; then
      filter_bam_by_min_aln_block "$raw_bam" "$bam" "$MIN_ALN_BLOCK"
      rm -f "$raw_bam"
    else
      mv "$raw_bam" "$bam"
    fi
    samtools index "$bam"
    mapped_reads=$(samtools view -c "${VIEW_FILTER_ARGS[@]}" "$bam" 2>/dev/null || echo 0)
    samtools depth "${DEPTH_ARGS[@]}" "$bam" > "$depthfile" || true
    read -r breadth meanD minD < <(awk -v L="$ref_len" -v mincov="$MIN_COV" '{d=$3; sum+=d; n+=1; if(d>=mincov) cov+=1; if(NR==1||d<min) min=d} END{if(L<=0||n==0){print "0 0 0"; exit} printf "%.4f %.2f %d\n", cov/L, sum/n, min}' "$depthfile")
    pct_mapped="0.00"
    if [[ "$reads_filtered" -gt 0 ]]; then pct_mapped="$(python3 - <<PY
rf=int("$reads_filtered"); mb=int("$mapped_reads")
print(f"{(100.0*mb/rf):.2f}")
PY
)"; fi
    echo -e "${sample_id}\t${sample_name}\t${barcode_folder}\t${bc}\t${reads_raw}\t${reads_trimmed}\t${reads_filtered}\t${ref_id}\t${ddns_group}\t${cluster}\t${ref_len}\t${mapped_reads}\t${pct_mapped}\t${breadth}\t${meanD}\t${minD}\t${bam}\t${depthfile}" >> "$sample_map_tsv"
  done < "$shortlist"

  cat <(tail -n +2 "$sample_map_tsv") >> "$BENCH_MASTER"

  candidate_file="$sample_dir/benchmark/candidate_refs.tsv"
  python3 - "$sample_map_tsv" "$candidate_file" "$RANK_MAP_TOPN" "$CONSENSUS_ALL_PASSING" "$ATTEMPT_CONSENSUS_ANY_MAPPED" "$MIN_BREADTH_FOR_CONS" "$MIN_MAPPED_FOR_CONS" "$MIN_MAPPED_PCENT_FOR_CONS" <<'PY'
import sys, csv
inp, outp, topn, allpassing, anymapped, minb, minm, minp = sys.argv[1:]
topn = int(topn); allpassing = int(allpassing); anymapped = int(anymapped); minb = float(minb); minm = int(minm); minp = float(minp)
rows=[]
with open(inp) as f:
    r=csv.DictReader(f, delimiter='	')
    for row in r:
        row['mapped_reads']=int(row['mapped_reads'])
        row['breadth']=float(row['breadth'])
        row['pct_mapped']=float(row['pct_mapped'])
        row['mean_depth']=float(row['mean_depth'])
        rows.append(row)
rows.sort(key=lambda x: (x['breadth'], x['pct_mapped'], x['mean_depth'], x['mapped_reads']), reverse=True)
passing=[r for r in rows if r['mapped_reads'] >= minm and r['breadth'] >= minb and r['pct_mapped'] >= minp]
if anymapped:
    chosen = [r for r in rows if r['mapped_reads'] > 0]
elif allpassing:
    chosen = passing
else:
    chosen = passing if topn == 0 else passing[:topn]
if not chosen and rows and rows[0]['mapped_reads'] > 0:
    chosen=[rows[0]]
with open(outp, 'w') as w:
    for r in chosen: w.write(r['ref_id'] + '\n')
PY

  while IFS= read -r cand_ref; do
    [[ -n "$cand_ref" ]] || continue
    meta_row="$(awk -F'\t' -v r="$cand_ref" 'NR>1 && $1==r {print; exit}' "$REF_META" || true)"
    [[ -n "$meta_row" ]] || continue
    ref_id="$cand_ref"
    safe_name="$(echo "$ref_id" | sed 's/[^A-Za-z0-9._-]/_/g')"
    ref_fa="$OUT/ref_cache/${safe_name}.fa"
    map_bam="$sample_dir/benchmark/${safe_name}.bam"
    medaka_bam="$sample_dir/consensus/${safe_name}.primary.bam"
    samtools view -b "${VIEW_FILTER_ARGS[@]}" "$map_bam" > "$medaka_bam" || true
    samtools index "$medaka_bam" || true
    hit_mapped=$(samtools view -c "${VIEW_FILTER_ARGS[@]}" "$medaka_bam" 2>/dev/null || echo 0)
    if [[ "$hit_mapped" -eq 0 ]]; then
      echo -e "${sample_id}\t${ref_id}\t1\tNO_MAPPED_READS\t\t0\t0\t0\t\t\t${medaka_bam}" >> "$sample_cons_tsv"
      continue
    fi
    cons_dir="$sample_dir/consensus/${safe_name}"
    mkdir -p "$cons_dir"
    if ! medaka_consensus -i "$medaka_bam" -d "$ref_fa" -o "$cons_dir" -t "$THREADS" -m "$RUN_MEDAKA_MODEL" > "$OUT/00_logs/${sample_id}.${safe_name}.medaka.log" 2>&1; then
      echo -e "${sample_id}\t${ref_id}\t1\tMEDAKA_FAIL\t\t0\t0\t0\t\t\t${medaka_bam}" >> "$sample_cons_tsv"
      continue
    fi
    cons_fa="$cons_dir/consensus.fasta"
    if [[ ! -s "$cons_fa" ]]; then
      echo -e "${sample_id}\t${ref_id}\t1\tNO_CONSENSUS\t\t0\t0\t0\t\t\t${medaka_bam}" >> "$sample_cons_tsv"
      continue
    fi
    lowbed="$cons_dir/lowcov.bed"
    depthfile="$sample_dir/benchmark/${safe_name}.depth.tsv"
    awk -v mincov="$MIN_COV" 'BEGIN{OFS="\t"} $3<mincov {s=$2-1; if(s<0)s=0; print $1,s,$2}' "$depthfile" > "$lowbed"
    masked="$cons_dir/consensus.masked.fasta"
    python3 - "$cons_fa" "$lowbed" > "$masked" <<'PY'
import sys
from collections import defaultdict
fa, bed = sys.argv[1], sys.argv[2]
seqs={}; name=None; buf=[]
with open(fa) as f:
    for line in f:
        line=line.rstrip('\n')
        if line.startswith('>'):
            if name is not None: seqs[name]=''.join(buf)
            name=line[1:].split()[0]; buf=[]
        else:
            buf.append(line)
    if name is not None: seqs[name]=''.join(buf)
mask=defaultdict(list)
with open(bed) as f:
    for line in f:
        if not line.strip(): continue
        chrom, st, en = line.rstrip('\n').split('\t')[:3]
        mask[chrom].append((int(st), int(en)))
def apply_mask(s, intervals):
    s=list(s)
    for st,en in intervals:
        st=max(st,0); en=min(en,len(s))
        for i in range(st,en): s[i]='N'
    return ''.join(s)
for chrom, seq in seqs.items():
    out = apply_mask(seq, mask.get(chrom, []))
    print(f'>{chrom}')
    for i in range(0, len(out), 80): print(out[i:i+80])
PY
    renamed="$cons_dir/consensus.masked.renamed.fasta"
    newhdr="${safe_sample}|${barcode_folder}|${ref_id}"
    awk -v h="$newhdr" 'BEGIN{done=0} /^>/{print ">"h; done=1; next} {print}' "$masked" > "$renamed"
    read -r cons_len n_count callable_bases muts mut_rate < <(python3 - "$renamed" "$ref_fa" <<'PY'
import sys
def read_seq(p):
    s=[]
    with open(p) as f:
        for line in f:
            if not line.startswith('>'): s.append(line.strip())
    return ''.join(s).upper()
c=read_seq(sys.argv[1]); r=read_seq(sys.argv[2]); n=min(len(c), len(r))
m=0; callable=0
for i in range(n):
    if c[i] == 'N': continue
    callable += 1
    if c[i] != r[i]: m += 1
ncount = c.count('N')
rate = '' if callable == 0 else f'{(m/callable):.6f}'
print(len(c), ncount, callable, m, rate)
PY
)
    echo -e "${sample_id}\t${ref_id}\t1\tOK\t${renamed}\t${cons_len}\t${n_count}\t${callable_bases}\t${muts}\t${mut_rate}\t${medaka_bam}" >> "$sample_cons_tsv"
  done < "$candidate_file"

  cat <(tail -n +2 "$sample_cons_tsv") >> "$CONS_MASTER"

  if [[ "$KEEP_ALL_BAMS" -eq 0 ]]; then
    find "$sample_dir/benchmark" -type f \( -name '*.bam' -o -name '*.bam.bai' \) -delete || true
    find "$sample_dir/screen" -type f \( -name '*.bam' -o -name '*.bam.bai' \) -delete || true
  fi
done

python3 - "$BENCH_MASTER" "$CONS_MASTER" "$OUT/benchmark_per_reference.tsv" "$OUT/summary.tsv" "$OUT/overall_reference_leaderboard.tsv" "$OUT/all_best_consensus.fasta" "$OUT/top_hits_report.tsv" "$FINAL_TOPK" "$MIN_CONS_LEN" "$MAX_N_FRAC" <<'PY'
import sys, csv, math, statistics, pathlib
map_tsv, cons_tsv, out_bench, out_sum, out_leader, out_best_cons, out_top_hits, final_topk, min_cons_len, max_n_frac = sys.argv[1:]
final_topk = int(final_topk)
min_cons_len = int(min_cons_len)
max_n_frac = float(max_n_frac)
def median_or_blank(vals):
    vals = [v for v in vals if v is not None]
    if not vals: return ''
    return statistics.median(vals)
cons = {}
with open(cons_tsv) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        key = (row['sample_id'], row['ref_id'])
        row['consensus_attempted'] = int(row['consensus_attempted'])
        row['consensus_len'] = int(row['consensus_len'] or 0)
        row['n_count'] = int(row['n_count'] or 0)
        row['callable_bases'] = int(row['callable_bases'] or 0)
        row['mutations'] = int(row['mutations'] or 0) if row['mutations'] else None
        row['mut_rate'] = float(row['mut_rate']) if row['mut_rate'] else None
        cons[key] = row
samples = {}; all_rows = []
with open(map_tsv) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        row['reads_raw'] = int(row['reads_raw']); row['reads_trimmed'] = int(row['reads_trimmed']); row['reads_filtered'] = int(row['reads_filtered'])
        row['ref_len'] = int(row['ref_len']); row['mapped_reads'] = int(row['mapped_reads']); row['pct_mapped'] = float(row['pct_mapped'])
        row['breadth'] = float(row['breadth']); row['mean_depth'] = float(row['mean_depth']); row['min_depth'] = int(float(row['min_depth']))
        c = cons.get((row['sample_id'], row['ref_id']))
        if c:
            row['consensus_attempted'] = c['consensus_attempted']; row['consensus_status_raw'] = c['consensus_status_raw']; row['consensus_fasta'] = c['consensus_fasta']
            row['consensus_len'] = c['consensus_len']; row['n_count'] = c['n_count']; row['callable_bases'] = c['callable_bases']
            row['mutations'] = c['mutations']; row['mut_rate'] = c['mut_rate']; row['medaka_bam'] = c['medaka_bam']
        else:
            row['consensus_attempted'] = 0; row['consensus_status_raw'] = 'NOT_ATTEMPTED'; row['consensus_fasta'] = ''
            row['consensus_len'] = 0; row['n_count'] = 0; row['callable_bases'] = 0; row['mutations'] = None; row['mut_rate'] = None; row['medaka_bam'] = ''
        row['n_frac'] = None
        if row['consensus_len'] > 0: row['n_frac'] = row['n_count'] / row['consensus_len']
        if row['consensus_attempted'] == 0: row['final_status'] = 'NO_CONSENSUS_ATTEMPT'
        elif row['consensus_status_raw'] != 'OK': row['final_status'] = row['consensus_status_raw']
        elif row['consensus_len'] < min_cons_len: row['final_status'] = 'CONSENSUS_SHORT'
        elif row['n_frac'] is not None and row['n_frac'] > max_n_frac: row['final_status'] = 'HIGH_N_FRAC'
        elif row['callable_bases'] == 0: row['final_status'] = 'NO_CALLABLE_BASES'
        else: row['final_status'] = 'CONSENSUS_PASS'
        samples.setdefault(row['sample_id'], []).append(row); all_rows.append(row)
for sample_id, rows in samples.items():
    rows_sorted_map = sorted(rows, key=lambda x: (x['breadth'], x['pct_mapped'], x['mean_depth'], x['mapped_reads']), reverse=True)
    for i, row in enumerate(rows_sorted_map, start=1): row['map_rank'] = i
    rows_sorted_final = sorted(rows, key=lambda x: (1 if x['final_status'] == 'CONSENSUS_PASS' else 0, x['breadth'], -(x['mut_rate'] if x['mut_rate'] is not None else math.inf), -(x['n_frac'] if x['n_frac'] is not None else math.inf), x['pct_mapped'], x['mean_depth'], x['mapped_reads']), reverse=True)
    for i, row in enumerate(rows_sorted_final, start=1):
        row['final_rank'] = i; row['is_best'] = 1 if i == 1 else 0; row['report_topk'] = 1 if i <= final_topk else 0
fieldnames = ['sample_id','sample_name','barcode_folder','bc','reads_raw','reads_trimmed','reads_filtered','ref_id','ddns_group','cluster','ref_len','mapped_reads','pct_mapped','breadth','mean_depth','min_depth','map_rank','consensus_attempted','consensus_status_raw','consensus_fasta','consensus_len','n_count','n_frac','callable_bases','mutations','mut_rate','final_status','final_rank','is_best','report_topk','map_bam','depth_tsv','medaka_bam']
with open(out_bench, 'w', newline='') as w:
    writer = csv.DictWriter(w, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    for sample_id in sorted(samples):
        for row in sorted(samples[sample_id], key=lambda x: x['final_rank']): writer.writerow(row)
top_hit_fields = ['sample_id','sample_name','barcode_folder','bc','hit_rank','ref_id','ddns_group','cluster','mapped_reads','breadth','mean_depth','min_depth','consensus_fasta','consensus_len','n_count','mutations','status']
with open(out_top_hits, 'w', newline='') as w:
    writer = csv.DictWriter(w, fieldnames=top_hit_fields, delimiter='	')
    writer.writeheader()
    for sample_id in sorted(samples):
        for row in sorted(samples[sample_id], key=lambda x: x['map_rank']):
            writer.writerow({'sample_id': row['sample_id'],'sample_name': row['sample_name'],'barcode_folder': row['barcode_folder'],'bc': row['bc'],'hit_rank': row['map_rank'],'ref_id': row['ref_id'],'ddns_group': row['ddns_group'],'cluster': row['cluster'],'mapped_reads': row['mapped_reads'],'breadth': f"{row['breadth']:.4f}",'mean_depth': f"{row['mean_depth']:.2f}",'min_depth': row['min_depth'],'consensus_fasta': row['consensus_fasta'],'consensus_len': row['consensus_len'],'n_count': row['n_count'],'mutations': '' if row['mutations'] is None else row['mutations'],'status': 'PASS' if row['final_status'] == 'CONSENSUS_PASS' else row['final_status']})
sum_fields = ['sample_id','sample_name','barcode_folder','bc','reads_raw','reads_trimmed','reads_filtered','best_ref','best_ddns_group','best_cluster','mapped_reads','pct_mapped','breadth','mean_depth','min_depth','consensus_len','n_count','n_frac','callable_bases','mutations','mut_rate','final_status','status','consensus_fasta']
best_rows = []
with open(out_sum, 'w', newline='') as w:
    writer = csv.DictWriter(w, fieldnames=sum_fields, delimiter='	'); writer.writeheader()
    for sample_id in sorted(samples):
        best = min(samples[sample_id], key=lambda x: x['final_rank'])
        best_rows.append(best)
        writer.writerow({'sample_id': best['sample_id'],'sample_name': best['sample_name'],'barcode_folder': best['barcode_folder'],'bc': best['bc'],'reads_raw': best['reads_raw'],'reads_trimmed': best['reads_trimmed'],'reads_filtered': best['reads_filtered'],'best_ref': best['ref_id'],'best_ddns_group': best['ddns_group'],'best_cluster': best['cluster'],'mapped_reads': best['mapped_reads'],'pct_mapped': f"{best['pct_mapped']:.2f}",'breadth': f"{best['breadth']:.4f}",'mean_depth': f"{best['mean_depth']:.2f}",'min_depth': best['min_depth'],'consensus_len': best['consensus_len'],'n_count': best['n_count'],'n_frac': '' if best['n_frac'] is None else f"{best['n_frac']:.6f}",'callable_bases': best['callable_bases'],'mutations': '' if best['mutations'] is None else best['mutations'],'mut_rate': '' if best['mut_rate'] is None else f"{best['mut_rate']:.6f}",'final_status': best['final_status'],'status': 'PASS' if best['final_status'] == 'CONSENSUS_PASS' else best['final_status'],'consensus_fasta': best['consensus_fasta']})
with open(out_best_cons, 'w') as w:
    for best in best_rows:
        p = best.get('consensus_fasta') or ''
        if p and pathlib.Path(p).exists(): w.write(pathlib.Path(p).read_text())
leader = {}
for row in all_rows:
    d = leader.setdefault(row['ref_id'], {'ddns_group': row['ddns_group'],'cluster': row['cluster'],'samples_seen': 0,'wins': 0,'consensus_passes': 0,'breadth': [],'pct_mapped': [],'mapped_reads': [],'mut_rate': [],'n_frac': [],'final_rank': []})
    d['samples_seen'] += 1; d['wins'] += 1 if row['is_best'] == 1 else 0; d['consensus_passes'] += 1 if row['final_status'] == 'CONSENSUS_PASS' else 0
    d['breadth'].append(row['breadth']); d['pct_mapped'].append(row['pct_mapped']); d['mapped_reads'].append(row['mapped_reads']); d['final_rank'].append(row['final_rank'])
    if row['mut_rate'] is not None: d['mut_rate'].append(row['mut_rate'])
    if row['n_frac'] is not None: d['n_frac'].append(row['n_frac'])
leader_fields = ['ref_id','ddns_group','cluster','samples_seen','wins','consensus_passes','median_final_rank','median_breadth','median_pct_mapped','median_mapped_reads','median_mut_rate','median_n_frac']
leader_rows = []
for ref_id, d in leader.items():
    leader_rows.append({'ref_id': ref_id,'ddns_group': d['ddns_group'],'cluster': d['cluster'],'samples_seen': d['samples_seen'],'wins': d['wins'],'consensus_passes': d['consensus_passes'],'median_final_rank': median_or_blank(d['final_rank']),'median_breadth': median_or_blank(d['breadth']),'median_pct_mapped': median_or_blank(d['pct_mapped']),'median_mapped_reads': median_or_blank(d['mapped_reads']),'median_mut_rate': median_or_blank(d['mut_rate']),'median_n_frac': median_or_blank(d['n_frac'])})
leader_rows.sort(key=lambda x: (x['wins'], x['consensus_passes'], -(x['median_final_rank'] if x['median_final_rank'] != '' else math.inf)), reverse=True)
with open(out_leader, 'w', newline='') as w:
    writer = csv.DictWriter(w, fieldnames=leader_fields, delimiter='\t')
    writer.writeheader()
    for row in leader_rows: writer.writerow(row)
PY

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$SCRIPT_DIR/report_tables.py" ]]; then
  python3 "$SCRIPT_DIR/report_tables.py" \
    --out "$OUT" \
    --fastq-pass "$FASTQ_PASS" \
    --barcode-map "$BARCODE_MAP" \
    --max-hits "$FINAL_TOPK" \
    --min-cons-len "$MIN_CONS_LEN" \
    --max-n-frac "$MAX_N_FRAC"
fi

if [[ -f "$SCRIPT_DIR/report_html.py" ]]; then
  extra=()
  if [[ -n "$RUN_REPORT_HTML" && -f "$RUN_REPORT_HTML" ]]; then
    extra+=( --run-report-html "$RUN_REPORT_HTML" )
  fi
  python3 "$SCRIPT_DIR/report_html.py" \
    --out "$OUT" \
    --run-name "$RUN_NAME" \
    "${extra[@]}"
fi

log "Creating results_bundle.zip"
python3 - "$OUT" <<'PY'
import sys, zipfile, pathlib
out = pathlib.Path(sys.argv[1])
zip_path = out / "results_bundle.zip"
include = []
for fn in [
    "summary.tsv",
    "benchmark_per_reference.tsv",
    "overall_reference_leaderboard.tsv",
    "all_best_consensus.fasta",
    "top_hits_report.tsv",
    "Table_Sample summary information.csv",
    "Table_Composition of samples.csv",
    "Table_Top hits.csv",
    "report.html",
    "00_logs/versions.txt",
]:
    p = out / fn
    if p.exists():
        include.append(p)
for p in (out / "per_sample").glob("*/consensus/*/consensus.masked.renamed.fasta"):
    include.append(p)
for p in (out / "per_sample").glob("*/idxstats.full.tsv"):
    include.append(p)
for p in (out / "00_logs").glob("*.log"):
    include.append(p)
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in include:
        if p.exists():
            z.write(p, arcname=str(p.relative_to(out)))
print(zip_path)
PY

log "Done"
log "Summary: $OUT/summary.tsv"
log "Top hits: $OUT/top_hits_report.tsv"
log "HTML report: $OUT/report.html"
log "Zip bundle: $OUT/results_bundle.zip"
