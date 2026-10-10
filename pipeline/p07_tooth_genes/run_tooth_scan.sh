#!/usr/bin/env bash
# run_tooth_scan.sh -- loi X9 / viec B5 (docs/05): quet MANH VO GENE RANG
# (ENAM, AMBN, AMELX, AMTN, DSPP, MMP20) va GULO trong bo gene chao mao bang
# TBLASTN voi protein ortholog ngoai nhom (ca sau My mac dinh, D06 docs/02)
# lam truy van; gene doi chung duong (MEPE, SPARC) di kem de kiem quy trinh.
#
# Thu tu (moi buoc resume theo "<output>.params co danh tinh input" nhu
# run_a2.sh / run_last.sh; --force chay lai het):
#   1 fetch-queries (python, host, CAN MANG NCBI)  -> data/b5/queries.faa + queries_manifest.tsv
#                                                    + 1 dong data/registry.csv (sha256)
#     (bo qua neu --queries FILE: dung FASTA truy van co san)
#   2 gunzip bo gene chao mao (gzip -dc tren host)  -> data/b5/genome.fna   (~1,2 GB, KHONG commit)
#   3 makeblastdb (docker|local)                     -> data/b5/blastdb/bulbul.*
#   4 tblastn (docker|local), -outfmt co dinh         -> data/b5/tblastn.tsv
#   5 parse (python)                                 -> loci.tsv, lesions.tsv, hits.tsv
#   6 extract (python, doc .fna.gz theo luong)        -> loci.fa, loci_translate.txt (kiem thu cong)
#   7 report (python)                                -> report.md (co canh bao bat buoc)
#   8 sha256sum cac dau ra nho                       -> SHA256SUMS.txt
#
# Chay tren Git Bash Windows (docker, nhu run_a2.sh/run_last.sh) HOAC WSL2/Linux
# (--mode local voi tblastn/makeblastdb tren PATH, vd conda env env/environment.yml).
#
# Docker image mac dinh: quay.io/biocontainers/blast:2.16.0--h66d330f_5
#   -> tag SUY TU bioconda repodata (goi blast-2.16.0-h66d330f_5; biocontainers
#      dat tag "<version>--<build>"), CHUA `docker pull` kiem trong lan viet
#      script nay (mang container bi chan). docs/02 ghi "can image BLAST, chua
#      duyet": chu du an duyet image roi pull; neu tag khong ton tai dung --image.
#
# KHONG co lenh xoa nao trong script nay (CLAUDE.md; su co 2026-09-08). Dau ra
# ghi de bang > / mv -f; thu muc cu khong bi don.
#
# CHUA chay tren du lieu that trong lan viet script (bo gene khong co o may
# viet code). Da qua: `bash -n`, smoke test --mode local tren bo gene tong hop
# nho (xem README.md muc "Smoke test").
set -euo pipefail

usage() {
  cat <<'EOF'
Dung: run_tooth_scan.sh [tuy chon]

Tien dieu kien: data/ncbi/GCA_013400435.1/..._genomic.fna.gz (D01, data/registry.csv).
Buoc 1 can mang toi NCBI va --email (hoac --queries FILE de bo buoc 1).

Tuy chon:
  --mode docker|local      docker (mac dinh, Git Bash + Docker Desktop) hoac
                           local (tblastn/makeblastdb tren PATH, WSL2/Linux).
  --image TAG              Docker image BLAST+ (xem dau file). Chi dung o mode docker.
  --email ADDR             Email cho NCBI E-utilities (bat buoc neu khong --queries).
  --api-key KEY            NCBI API key (tuy chon).
  --organism "Ten loai"    Ngoai nhom lay protein truy van (mac dinh "Alligator mississippiensis").
  --genes A,B,...          Gene muc tieu (mac dinh ENAM,AMBN,AMELX,AMTN,DSPP,MMP20,GULO).
  --controls A,B,...       Gene doi chung duong (mac dinh MEPE,SPARC).
  --queries FILE           FASTA protein truy van co san (bo buoc fetch; neu header
                           khong theo dang GENE|ROLE|ACC thi gene = id, role = target).
  --manifest FILE          queries_manifest.tsv di kem --queries (tuy chon).
  --genome FILE            Bo gene .fna.gz (mac dinh D01 trong --data-dir).
  --evalue E               Nguong e-value tblastn + parse (mac dinh 1e-5).
  --threads N              -num_threads cho tblastn (mac dinh so CPU - 2, hoac 4).
  --pad N                  So nt them moi ben khi trich locus (mac dinh 2000).
  --force                  Chay lai moi buoc du output + tham so da khop.
  --data-dir DIR           Thu muc data/ (mac dinh <repo>/data).
  -h, --help               In huong dan nay va thoat (khong dung docker/mang).
EOF
}

MODE="${MODE:-docker}"
IMG_BLAST="${IMG_BLAST:-quay.io/biocontainers/blast:2.16.0--h66d330f_5}"
EMAIL="${NCBI_EMAIL:-}"
API_KEY="${NCBI_API_KEY:-}"
ORGANISM="${ORGANISM:-Alligator mississippiensis}"
GENES="${GENES:-ENAM,AMBN,AMELX,AMTN,DSPP,MMP20,GULO}"
CONTROLS="${CONTROLS:-MEPE,SPARC}"
QUERIES_OVERRIDE=""
MANIFEST_OVERRIDE=""
GENOME_OVERRIDE=""
EVALUE="${EVALUE:-1e-5}"
THREADS="${THREADS:-}"
PAD="${PAD:-2000}"
FORCE=0
DATA_DIR_OVERRIDE=""

while [ $# -gt 0 ]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --image) IMG_BLAST="$2"; shift 2 ;;
    --email) EMAIL="$2"; shift 2 ;;
    --api-key) API_KEY="$2"; shift 2 ;;
    --organism) ORGANISM="$2"; shift 2 ;;
    --genes) GENES="$2"; shift 2 ;;
    --controls) CONTROLS="$2"; shift 2 ;;
    --queries) QUERIES_OVERRIDE="$2"; shift 2 ;;
    --manifest) MANIFEST_OVERRIDE="$2"; shift 2 ;;
    --genome) GENOME_OVERRIDE="$2"; shift 2 ;;
    --evalue) EVALUE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --pad) PAD="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --data-dir) DATA_DIR_OVERRIDE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "run_tooth_scan.sh: tuy chon khong biet: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$MODE" in docker|local) ;; *) echo "run_tooth_scan.sh: --mode phai la docker hoac local" >&2; exit 2 ;; esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Script nam sau 2 cap (pipeline/p07_tooth_genes/) -> len 2 cap. Lesson 2026-09-09:
# run_last.sh tung dung "$SCRIPT_DIR/.." va tro sai vao pipeline/data.
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="${DATA_DIR_OVERRIDE:-$PROJECT_DIR/data}"
OUT_DIR="$DATA_DIR/b5"
LOG="$OUT_DIR/run.log"
PY="$SCRIPT_DIR/tooth_scan.py"
REGISTRY="$DATA_DIR/registry.csv"
GENOME_GZ="${GENOME_OVERRIDE:-$DATA_DIR/ncbi/GCA_013400435.1/GCA_013400435.1_ASM1340043v1_genomic.fna.gz}"

if command -v python3 >/dev/null 2>&1; then PY_BIN=python3; elif command -v python >/dev/null 2>&1; then PY_BIN=python; else
  echo "run_tooth_scan.sh: khong tim thay python3/python" >&2; exit 3; fi

mkdir -p "$OUT_DIR"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

# --- Ham dung chung, sao chep NGUYEN VAN tu p02_core_map/run_a2.sh (qua
# p02b_last/run_last.sh; khong import cheo giua cac thu muc pipeline). ---------

params_file_for() {
  # params_file_for FILE -> ten file ".params" canh FILE.
  echo "${1}.params"
}

input_id() {
  # input_id FILE... -> chuoi "ten:size:mtime" cho tung input, de PARAMS_STR
  # doi khi INPUT doi (khong chi khi tham so doi).
  local out="" f
  for f in "$@"; do
    if [ -f "$f" ]; then
      out="$out $(basename "$f"):$(stat -c %s "$f" 2>/dev/null):$(stat -c %Y "$f" 2>/dev/null)"
    else
      out="$out $(basename "$f"):missing"
    fi
  done
  echo "${out# }"
}

step_needed() {
  # step_needed PARAMS_FILE PARAMS_STR FILE... -> "can chay" (return 0) neu
  # FORCE=1, thieu 1 FILE bat ky, HOAC PARAMS_STR khac lan chay truoc.
  local pf="$1" ps="$2"
  shift 2
  if [ "$FORCE" -eq 1 ]; then
    return 0
  fi
  for f in "$@"; do
    [ -f "$f" ] || return 0
  done
  if [ -f "$pf" ] && [ "$(cat "$pf")" = "$ps" ]; then
    return 1
  fi
  return 0
}

save_params() {
  printf '%s' "$2" > "$1"
}

# --- Docker / local ---------------------------------------------------------

win_path() {
  # Git Bash: doi sang duong dan Windows cho -v (MSYS_NO_PATHCONV=1 tat tu quy
  # doi). Linux/WSL2 co docker: khong co cygpath -> giu nguyen duong dan.
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -w "$1"
  else
    echo "$1"
  fi
}

CPU_COUNT="$(nproc 2>/dev/null || echo 0)"
if [ "$CPU_COUNT" -gt 2 ] 2>/dev/null; then DEFAULT_THREADS=$((CPU_COUNT - 2)); else DEFAULT_THREADS=4; fi
THREADS="${THREADS:-$DEFAULT_THREADS}"

if [ "$MODE" = "docker" ]; then
  MOUNT_DATA_DIR="$(win_path "$DATA_DIR")"
  blast_cmd() {
    # blast_cmd <tool> <args...>: duong dan trong container bat dau bang /data
    # (= DATA_DIR). Moi duong dan file duoi day deu viet theo /data/b5/... nen
    # --data-dir doi tuy y van dung (lesson 2026-09-09 ve duong dan cung).
    MSYS_NO_PATHCONV=1 docker run --rm -v "${MOUNT_DATA_DIR}:/data" "$IMG_BLAST" "$@"
  }
  C_GENOME="/data/b5/genome.fna"
  C_DB="/data/b5/blastdb/bulbul"
  C_QUERIES="/data/b5/queries.faa"
else
  for t in makeblastdb tblastn; do
    command -v "$t" >/dev/null 2>&1 || { echo "run_tooth_scan.sh: --mode local can '$t' tren PATH (conda install blast, xem env/environment.yml)" >&2; exit 3; }
  done
  blast_cmd() { "$@"; }
  C_GENOME="$OUT_DIR/genome.fna"
  C_DB="$OUT_DIR/blastdb/bulbul"
  C_QUERIES="$OUT_DIR/queries.faa"
fi

if [ ! -f "$GENOME_GZ" ]; then
  echo "run_tooth_scan.sh: thieu bo gene $GENOME_GZ (D01 trong data/registry.csv) - tai truoc hoac dung --genome" >&2
  exit 4
fi

QUERIES="$OUT_DIR/queries.faa"
MANIFEST="$OUT_DIR/queries_manifest.tsv"
GENOME_FNA="$OUT_DIR/genome.fna"
DB_DIR="$OUT_DIR/blastdb"
DB_MARKER="$DB_DIR/bulbul.nsq"
TBLASTN_OUT="$OUT_DIR/tblastn.tsv"
LOCI="$OUT_DIR/loci.tsv"
LESIONS="$OUT_DIR/lesions.tsv"
HITS="$OUT_DIR/hits.tsv"
LOCI_FA="$OUT_DIR/loci.fa"
LOCI_TR="$OUT_DIR/loci_translate.txt"
REPORT="$OUT_DIR/report.md"
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sframe qseq sseq"

log "=== run_tooth_scan.sh bat dau (mode=$MODE force=$FORCE threads=$THREADS evalue=$EVALUE organism=\"$ORGANISM\" genes=$GENES controls=$CONTROLS queries_override=${QUERIES_OVERRIDE:-none}) ==="

# --- Buoc 1: truy van protein ngoai nhom ------------------------------------
if [ -n "$QUERIES_OVERRIDE" ]; then
  if [ ! -f "$QUERIES_OVERRIDE" ]; then echo "run_tooth_scan.sh: --queries $QUERIES_OVERRIDE khong ton tai" >&2; exit 4; fi
  Q_PARAMS_FILE="$(params_file_for "$QUERIES")"
  Q_PARAMS_STR="copy_from=$(input_id "$QUERIES_OVERRIDE")"
  if step_needed "$Q_PARAMS_FILE" "$Q_PARAMS_STR" "$QUERIES"; then
    log "STEP queries_copy START ($QUERIES_OVERRIDE -> $QUERIES)"
    cp -f "$QUERIES_OVERRIDE" "$QUERIES"
    if [ -n "$MANIFEST_OVERRIDE" ]; then cp -f "$MANIFEST_OVERRIDE" "$MANIFEST"; fi
    save_params "$Q_PARAMS_FILE" "$Q_PARAMS_STR"
    log "STEP queries_copy END"
  else
    log "STEP queries_copy SKIP (da khop; --force de chay lai)"
  fi
else
  if [ -z "$EMAIL" ]; then echo "run_tooth_scan.sh: can --email (NCBI E-utilities) hoac --queries FILE" >&2; exit 2; fi
  Q_PARAMS_FILE="$(params_file_for "$QUERIES")"
  Q_PARAMS_STR="organism=$ORGANISM genes=$GENES controls=$CONTROLS"
  if step_needed "$Q_PARAMS_FILE" "$Q_PARAMS_STR" "$QUERIES" "$MANIFEST"; then
    log "STEP fetch_queries START"
    API_ARGS=()
    if [ -n "$API_KEY" ]; then API_ARGS=(--api-key "$API_KEY"); fi
    REG_ARGS=()
    if [ -f "$REGISTRY" ]; then REG_ARGS=(--registry "$REGISTRY"); else log "canh bao: khong thay $REGISTRY, khong ghi registry"; fi
    "$PY_BIN" -B "$PY" fetch-queries --genes "$GENES" --controls "$CONTROLS" --organism "$ORGANISM" \
      --email "$EMAIL" "${API_ARGS[@]}" --out "$QUERIES" --manifest "$MANIFEST" "${REG_ARGS[@]}" 2>&1 | tee -a "$LOG"
    save_params "$Q_PARAMS_FILE" "$Q_PARAMS_STR"
    log "STEP fetch_queries END"
  else
    log "STEP fetch_queries SKIP (da khop; --force de chay lai)"
  fi
fi
if [ ! -s "$QUERIES" ]; then echo "run_tooth_scan.sh: $QUERIES rong - khong co protein truy van nao" >&2; exit 5; fi

# --- Buoc 2: giai nen bo gene (makeblastdb khong doc .gz) -------------------
G_PARAMS_FILE="$(params_file_for "$GENOME_FNA")"
G_PARAMS_STR="in=[$(input_id "$GENOME_GZ")]"
if step_needed "$G_PARAMS_FILE" "$G_PARAMS_STR" "$GENOME_FNA"; then
  log "STEP gunzip_genome START ($GENOME_GZ -> $GENOME_FNA)"
  gzip -dc "$GENOME_GZ" > "$GENOME_FNA.tmp"
  mv -f "$GENOME_FNA.tmp" "$GENOME_FNA"
  save_params "$G_PARAMS_FILE" "$G_PARAMS_STR"
  log "STEP gunzip_genome END size=$(stat -c %s "$GENOME_FNA")"
else
  log "STEP gunzip_genome SKIP (da khop; --force de chay lai)"
fi

# --- Buoc 3: makeblastdb ----------------------------------------------------
mkdir -p "$DB_DIR"
D_PARAMS_FILE="$(params_file_for "$DB_MARKER")"
D_PARAMS_STR="mode=$MODE image=$IMG_BLAST in=[$(input_id "$GENOME_FNA")]"
if step_needed "$D_PARAMS_FILE" "$D_PARAMS_STR" "$DB_MARKER"; then
  log "STEP makeblastdb START"
  blast_cmd makeblastdb -in "$C_GENOME" -dbtype nucl -out "$C_DB" -title bulbul_GCA_013400435.1 2>&1 | tee -a "$LOG"
  save_params "$D_PARAMS_FILE" "$D_PARAMS_STR"
  log "STEP makeblastdb END"
else
  log "STEP makeblastdb SKIP (da khop; --force de chay lai)"
fi

# --- Buoc 4: tblastn --------------------------------------------------------
T_PARAMS_FILE="$(params_file_for "$TBLASTN_OUT")"
T_PARAMS_STR="mode=$MODE image=$IMG_BLAST evalue=$EVALUE in=[$(input_id "$QUERIES" "$DB_MARKER")]"
if step_needed "$T_PARAMS_FILE" "$T_PARAMS_STR" "$TBLASTN_OUT"; then
  log "STEP tblastn START (threads=$THREADS evalue=$EVALUE)"
  t0=$(date +%s)
  # -seg no: khong che vung do phuc tap thap (protein men rang giau P/Q/S rat
  # lap -> seg mac dinh se che mat chinh vung can tim). -max_target_seqs 50:
  # gia gene thuong vo thanh nhieu scaffold o lap rap muc scaffold.
  blast_cmd tblastn -query "$C_QUERIES" -db "$C_DB" -evalue "$EVALUE" -seg no \
    -max_target_seqs 50 -num_threads "$THREADS" -outfmt "$OUTFMT" > "$TBLASTN_OUT.tmp"
  mv -f "$TBLASTN_OUT.tmp" "$TBLASTN_OUT"
  save_params "$T_PARAMS_FILE" "$T_PARAMS_STR"
  log "STEP tblastn END elapsed=$(( $(date +%s) - t0 ))s hsps=$(wc -l < "$TBLASTN_OUT")"
else
  log "STEP tblastn SKIP (da khop; --force de chay lai)"
fi

# --- Buoc 5-7: parse, extract, report (python) -------------------------------
MAN_ARGS=()
if [ -f "$MANIFEST" ]; then MAN_ARGS=(--manifest "$MANIFEST"); fi

P_PARAMS_FILE="$(params_file_for "$LOCI")"
P_PARAMS_STR="evalue=$EVALUE in=[$(input_id "$TBLASTN_OUT" "$MANIFEST")]"
if step_needed "$P_PARAMS_FILE" "$P_PARAMS_STR" "$LOCI" "$LESIONS" "$HITS"; then
  log "STEP parse START"
  "$PY_BIN" -B "$PY" parse --tblastn "$TBLASTN_OUT" "${MAN_ARGS[@]}" --max-evalue "$EVALUE" \
    --out-loci "$LOCI" --out-lesions "$LESIONS" --out-hits "$HITS" 2>&1 | tee -a "$LOG"
  save_params "$P_PARAMS_FILE" "$P_PARAMS_STR"
  log "STEP parse END"
else
  log "STEP parse SKIP (da khop; --force de chay lai)"
fi

E_PARAMS_FILE="$(params_file_for "$LOCI_FA")"
E_PARAMS_STR="pad=$PAD in=[$(input_id "$LOCI" "$GENOME_GZ")]"
if step_needed "$E_PARAMS_FILE" "$E_PARAMS_STR" "$LOCI_FA" "$LOCI_TR"; then
  log "STEP extract START (pad=$PAD)"
  if [ "$(wc -l < "$LOCI")" -gt 1 ]; then
    "$PY_BIN" -B "$PY" extract --fna "$GENOME_GZ" --loci "$LOCI" --pad "$PAD" --out "$LOCI_FA" --translate "$LOCI_TR" 2>&1 | tee -a "$LOG"
  else
    log "extract: loci.tsv khong co locus nao -> ghi file rong"
    : > "$LOCI_FA"; : > "$LOCI_TR"
  fi
  save_params "$E_PARAMS_FILE" "$E_PARAMS_STR"
  log "STEP extract END"
else
  log "STEP extract SKIP (da khop; --force de chay lai)"
fi

R_PARAMS_FILE="$(params_file_for "$REPORT")"
R_PARAMS_STR="in=[$(input_id "$LOCI" "$LESIONS" "$MANIFEST")]"
if step_needed "$R_PARAMS_FILE" "$R_PARAMS_STR" "$REPORT"; then
  log "STEP report START"
  "$PY_BIN" -B "$PY" report --loci "$LOCI" --lesions "$LESIONS" "${MAN_ARGS[@]}" --out "$REPORT" 2>&1 | tee -a "$LOG"
  save_params "$R_PARAMS_FILE" "$R_PARAMS_STR"
  log "STEP report END"
else
  log "STEP report SKIP (da khop; --force de chay lai)"
fi

# --- Buoc 8: sha256 cac dau ra nho (khong tinh genome.fna / blastdb) ----------
log "STEP sha256 START"
( cd "$OUT_DIR" && sha256sum queries.faa tblastn.tsv loci.tsv lesions.tsv hits.tsv loci.fa loci_translate.txt report.md $( [ -f queries_manifest.tsv ] && echo queries_manifest.tsv ) > SHA256SUMS.txt )
log "STEP sha256 END"

log "=== run_tooth_scan.sh xong. Doc $REPORT (muc cuoi truoc), roi kiem thu cong $LOCI_TR ==="
