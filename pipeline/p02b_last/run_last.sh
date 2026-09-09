#!/usr/bin/env bash
# run_last.sh -- loi X4: phan tu bao ton NGAN (<100 bp) khong duoc minimap2
# can (loi A2, xem ../p02_core_map/run_a2.sh) -> can lai bang LAST theo tung
# manh bo gene chao mao -> BED toa do chao mao -> gop voi bulbul_core.bed
# (dau ra paf2bed cua A2) thanh bulbul_core_all.bed.
#
# Chay tren Git Bash Windows. Thu tu:
#   filter (python)                          elements.fa+elements.bed -> short.fa
#   split_fasta.py (python, tai dung tu p02_core_map)  -> genome_chunks_last_<N>/chunk_NN.fa
#   [moi manh] docker: lastdb -uNEAR  ->  docker: lastal -D1e6 | docker: last-split
#     -> noi MAF cua moi manh lai thanh aln.maf
#   maf2bed (python)                         aln.maf -> bulbul_short.bed
#   merge-beds (python)                      bulbul_core.bed (A2) + bulbul_short.bed -> bulbul_core_all.bed
#
# Vi sao can buoc nay: minimap2 (ke ca tham so seed da toi uu cho A2, xem
# run_a2.sh) chi can duoc 0,3% phan tu 20-49bp va 7,8% phan tu 50-99bp tu ga
# sang chao mao (do that, docs/04 muc 10b) - nhom nay chiem ~80% so phan tu
# va ~41% tong bp bao ton sau merge. LAST la cong cu chuan cho can trinh tu
# ngan giua loai xa nhau. Xem research/briefs/X4-last-short-elements.md.
#
# Moi buoc ghi START/END/SKIP kem gio vao data/a2_last/run.log, dung cung co
# che "<output>.params co danh tinh input" nhu run_a2.sh (ham input_id/
# step_needed/save_params sao chep NGUYEN VAN tu run_a2.sh - xem README.md
# cung thu muc: khong sua run_a2.sh, khong import cheo giua 2 thu muc pipeline,
# dung quy uoc "moi script tu chua" san co cua repo). --force bo qua ca 2
# dieu kien, luon chay lai.
#
# KHAC run_a2.sh o 1 diem CO CHU DICH: buoc chia manh KHONG xoa thu muc manh
# cu truoc khi chia lai (run_a2.sh dung `rm -rf` cho viec nay). Du an cam
# tuyet doi lenh xoa trong script (xem CLAUDE.md + tasks/RESTORE_LOG-2026-09-08.txt
# - su co xoa nham 2026-09-08). Thay vao do, thu muc manh duoc dat ten KEM
# --chunk-bases (genome_chunks_last_<N>) nen doi --chunk-bases se tu chuyen
# sang thu muc moi (rong) thay vi de lai manh cu cung ten voi so luong khac -
# khong con nguy co manh "mo côi" tu lan chia truoc voi --chunk-bases khac.
# Doi lai: doi --chunk-bases nhieu lan se giu ca cac thu muc manh cu tren
# dia (khong tu don) - xem README.md.
#
# CANH BAO BO NHO (doc truoc khi chay tren VM Docker ~2GB RAM): lastdb voi
# -uNEAR cho 1 manh 90 Mbp (mac dinh --chunk-bases) co the vuot 2GB va bi
# container giet (exit 137). Neu gap loi nay, chay lai voi --chunk-bases
# 45000000.
#
# Docker image (tag dung nguyen van, da pull san - xem research/briefs/X4-*.md):
#   quay.io/biocontainers/last:1654--h5814d7d_1  (lastdb, lastal, last-split)
#
# CHUA chay docker/LAST that trong lan viet script nay (theo dung pham vi
# duoc giao) - script moi chi qua kiem tra `bash -n` (cu phap), chua qua
# chay thuc te tren du lieu that.
set -euo pipefail

usage() {
  cat <<'EOF'
Dung: run_last.sh [tuy chon]

Tien dieu kien: da chay p02_core_map/run_a2.sh toi it nhat buoc twoBitToFa,
tuc la co san data/a2/elements.bed va data/a2/elements.fa. Buoc merge-beds
cuoi cung con can data/a2/bulbul_core.bed (dau ra paf2bed cua A2).

Tuy chon:
  --force                  Chay lai moi buoc du output (va tham so) da khop.
  --threads N               So luong thread cho lastdb/lastal (mac dinh: so
                           CPU - 2, hoac 12 neu khong doc duoc so CPU).
  --chunk-bases N            So base toi da moi manh khi chia rieng bo gene
                           chao mao cho LAST (mac dinh 90000000 = 90 Mbp -
                           NHO HON mac dinh 180000000 cua run_a2.sh vi lastdb
                           -uNEAR ton RAM hon minimap2; giam xuong 45000000
                           neu container bi giet - exit 137). Moi gia tri
                           dung 1 thu muc manh rieng (xem ghi chu tren dau file).
  --max-len N                Do dai toi da (bp) de 1 phan tu duoc coi la
                           "ngan" va dua vao LAST (mac dinh 99; lenh con
                           filter cua short_elements.py giu --min-len mac
                           dinh rieng = 20, khong cau hinh o day).
  --lastdb-extra "CHUOI"    Tham so them SAU `-uNEAR` cua moi lenh lastdb
                           (mac dinh rong).
  --lastal-extra "CHUOI"    Tham so them SAU `-D1e6` cua moi lenh lastal
                           (mac dinh rong).
  --data-dir DIR            Thu muc data/ (mac dinh <repo>/data) - giong
                           run_a2.sh, dung khi test/doi vi tri du lieu.
  -h, --help                In huong dan nay va thoat (khong dung docker).
EOF
}

FORCE=0
THREADS="${THREADS:-}"
CHUNK_BASES="${CHUNK_BASES:-90000000}"
MAX_LEN="${MAX_LEN:-99}"
LASTDB_EXTRA="${LASTDB_EXTRA:-}"
LASTAL_EXTRA="${LASTAL_EXTRA:-}"
DATA_DIR_OVERRIDE=""

while [ $# -gt 0 ]; do
  case "$1" in
    --force) FORCE=1; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --chunk-bases) CHUNK_BASES="$2"; shift 2 ;;
    --max-len) MAX_LEN="$2"; shift 2 ;;
    --lastdb-extra) LASTDB_EXTRA="$2"; shift 2 ;;
    --lastal-extra) LASTAL_EXTRA="$2"; shift 2 ;;
    --data-dir) DATA_DIR_OVERRIDE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "run_last.sh: tuy chon khong biet: $1" >&2; usage >&2; exit 2 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="${DATA_DIR_OVERRIDE:-$PROJECT_DIR/data}"
CORE_MAP_DIR="$SCRIPT_DIR/../p02_core_map"
A2_DIR="$DATA_DIR/a2"
OUT_DIR="$DATA_DIR/a2_last"
LOG="$OUT_DIR/run.log"
PY="$SCRIPT_DIR/short_elements.py"
SPLIT_PY="$CORE_MAP_DIR/split_fasta.py"

ELEMENTS_BED="$A2_DIR/elements.bed"
ELEMENTS_FA="$A2_DIR/elements.fa"
BULBUL_CORE_BED="$A2_DIR/bulbul_core.bed"
BULBUL_FNA_HOST="$DATA_DIR/ncbi/GCA_013400435.1/GCA_013400435.1_ASM1340043v1_genomic.fna.gz"

IMG_LAST="quay.io/biocontainers/last:1654--h5814d7d_1"

mkdir -p "$OUT_DIR"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"
}

# --- Ham dung chung, sao chep NGUYEN VAN tu p02_core_map/run_a2.sh (khong
# sua run_a2.sh, khong import cheo - xem ghi chu dau file). ------------------

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

win_path() {
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -w "$1"
  else
    echo "run_last.sh: khong tim thay cygpath (can Git Bash / MSYS2) de doi duong dan cho docker -v" >&2
    exit 3
  fi
}

WIN_DATA_DIR="$(win_path "$DATA_DIR")"

CPU_COUNT="$(nproc 2>/dev/null || echo 0)"
if [ "$CPU_COUNT" -gt 2 ] 2>/dev/null; then
  DEFAULT_THREADS=$((CPU_COUNT - 2))
else
  DEFAULT_THREADS=12
fi
THREADS="${THREADS:-$DEFAULT_THREADS}"

if [ ! -f "$ELEMENTS_BED" ] || [ ! -f "$ELEMENTS_FA" ]; then
  echo "run_last.sh: thieu $ELEMENTS_BED hoac $ELEMENTS_FA - chay p02_core_map/run_a2.sh truoc (can toi buoc twoBitToFa)" >&2
  exit 4
fi

SHORT_FA="$OUT_DIR/short.fa"
FILTER_STATS="$OUT_DIR/filter_stats.tsv"
GENOME_CHUNKS_DIR="$OUT_DIR/genome_chunks_last_${CHUNK_BASES}"
SPLIT_MARKER="$GENOME_CHUNKS_DIR/.split_done"
LASTDB_DIR="$OUT_DIR/lastdb"
MAF_OUT="$OUT_DIR/aln.maf"
BULBUL_SHORT_BED="$OUT_DIR/bulbul_short.bed"
UNMAPPED_SHORT="$OUT_DIR/unmapped_short.txt"
MAP_STATS_SHORT="$OUT_DIR/map_stats_short.tsv"
BULBUL_CORE_ALL_BED="$OUT_DIR/bulbul_core_all.bed"
MERGE_BEDS_STATS="$OUT_DIR/merge_beds_stats.tsv"

log "=== run_last.sh bat dau (force=$FORCE threads=$THREADS chunk-bases=$CHUNK_BASES max-len=$MAX_LEN lastdb-extra=\"$LASTDB_EXTRA\" lastal-extra=\"$LASTAL_EXTRA\") ==="

# --- Buoc 1: filter (loc phan tu ngan ra short.fa) --------------------------
FILTER_PARAMS_FILE="$(params_file_for "$SHORT_FA")"
FILTER_PARAMS_STR="max_len=$MAX_LEN in=[$(input_id "$ELEMENTS_FA" "$ELEMENTS_BED")]"
if step_needed "$FILTER_PARAMS_FILE" "$FILTER_PARAMS_STR" "$SHORT_FA" "$FILTER_STATS"; then
  log "STEP filter START"
  python -B "$PY" filter --in-fa "$ELEMENTS_FA" --in-bed "$ELEMENTS_BED" \
    --out "$SHORT_FA" --max-len "$MAX_LEN" --stats "$FILTER_STATS"
  save_params "$FILTER_PARAMS_FILE" "$FILTER_PARAMS_STR"
  log "STEP filter END"
else
  log "STEP filter SKIP (output + tham so da khop; dung --force de chay lai)"
fi

# --- Buoc 2: chia bo gene chao mao RIENG cho LAST (thu muc theo chunk-bases,
# xem ghi chu dau file - KHONG xoa thu muc manh cu). -------------------------
SPLIT_PARAMS_FILE="$(params_file_for "$SPLIT_MARKER")"
SPLIT_PARAMS_STR="chunk_bases=$CHUNK_BASES in=[$(input_id "$BULBUL_FNA_HOST")]"
if step_needed "$SPLIT_PARAMS_FILE" "$SPLIT_PARAMS_STR" "$SPLIT_MARKER"; then
  log "STEP split_genome_last START (chunk-bases=$CHUNK_BASES, thu muc $GENOME_CHUNKS_DIR)"
  mkdir -p "$GENOME_CHUNKS_DIR"
  python -B "$SPLIT_PY" --in "$BULBUL_FNA_HOST" --out-dir "$GENOME_CHUNKS_DIR" --max-bases "$CHUNK_BASES" 2>&1 | tee -a "$LOG"
  touch "$SPLIT_MARKER"
  save_params "$SPLIT_PARAMS_FILE" "$SPLIT_PARAMS_STR"
  log "STEP split_genome_last END"
else
  log "STEP split_genome_last SKIP (da co manh voi cung chunk-bases; dung --force de chia lai - ghi de cung thu muc, khong xoa)"
fi

# --- Buoc 3: lastdb + lastal | last-split tung manh, noi MAF ----------------
LAST_PARAMS_FILE="$(params_file_for "$MAF_OUT")"
LAST_PARAMS_STR="chunk_bases=$CHUNK_BASES lastdb_extra=$LASTDB_EXTRA lastal_extra=$LASTAL_EXTRA threads=$THREADS in=[$(input_id "$SHORT_FA" "$BULBUL_FNA_HOST")]"
if step_needed "$LAST_PARAMS_FILE" "$LAST_PARAMS_STR" "$MAF_OUT"; then
  log "STEP last_chunks START (threads=$THREADS lastdb-extra=\"$LASTDB_EXTRA\" lastal-extra=\"$LASTAL_EXTRA\")"
  MAF_TMP="$MAF_OUT.tmp"
  : > "$MAF_TMP"
  mkdir -p "$LASTDB_DIR"
  shopt -s nullglob
  chunks=("$GENOME_CHUNKS_DIR"/chunk_*.fa)
  shopt -u nullglob
  if [ "${#chunks[@]}" -eq 0 ]; then
    echo "run_last.sh: khong tim thay manh nao trong $GENOME_CHUNKS_DIR (split_genome_last that bai?)" >&2
    exit 5
  fi
  for chunk in "${chunks[@]}"; do
    cb="$(basename "$chunk" .fa)"
    t0=$(date +%s)
    log "CHUNK $cb START"
    mkdir -p "$LASTDB_DIR/$cb"
    # lastdb: xay chi muc LAST cho 1 manh bo gene chao mao (-uNEAR: seed
    # phu hop trinh tu gan giong nhau - xem brief X4 muc "Viec" #2).
    MSYS_NO_PATHCONV=1 docker run --rm -v "${WIN_DATA_DIR}:/data" "$IMG_LAST" \
        lastdb -P "$THREADS" -uNEAR $LASTDB_EXTRA \
        "/data/a2_last/lastdb/$cb/db" "/data/a2_last/genome_chunks_last_${CHUNK_BASES}/$cb.fa"
    # lastal | last-split: can short.fa vao chi muc vua xay, roi giai quyet
    # cac hit chong lan/split (2 container rieng noi qua pipe tren host vi
    # last-split doc MAF tu stdin, khong nhan tham so file dau vao).
    MSYS_NO_PATHCONV=1 docker run --rm -v "${WIN_DATA_DIR}:/data" "$IMG_LAST" \
        lastal -P "$THREADS" -D1e6 $LASTAL_EXTRA \
        "/data/a2_last/lastdb/$cb/db" "/data/a2_last/short.fa" \
      | MSYS_NO_PATHCONV=1 docker run --rm -i -v "${WIN_DATA_DIR}:/data" "$IMG_LAST" \
        last-split \
      >> "$MAF_TMP"
    t1=$(date +%s)
    log "CHUNK $cb END elapsed=$((t1 - t0))s lines=$(wc -l < "$MAF_TMP")"
  done
  mv -f "$MAF_TMP" "$MAF_OUT"
  save_params "$LAST_PARAMS_FILE" "$LAST_PARAMS_STR"
  log "STEP last_chunks END lines=$(wc -l < "$MAF_OUT")"
else
  log "STEP last_chunks SKIP (output + tham so da khop; dung --force de chay lai)"
fi

# --- Buoc 4: maf2bed --------------------------------------------------------
MAF2BED_PARAMS_FILE="$(params_file_for "$BULBUL_SHORT_BED")"
MAF2BED_PARAMS_STR="in=[$(input_id "$MAF_OUT" "$ELEMENTS_BED")]"
if step_needed "$MAF2BED_PARAMS_FILE" "$MAF2BED_PARAMS_STR" "$BULBUL_SHORT_BED" "$UNMAPPED_SHORT" "$MAP_STATS_SHORT"; then
  log "STEP maf2bed START"
  python -B "$PY" maf2bed --maf "$MAF_OUT" --out "$BULBUL_SHORT_BED" \
    --unmapped "$UNMAPPED_SHORT" --stats "$MAP_STATS_SHORT" \
    --min-identity 0.70 --min-coverage 0.80 --query-list "$ELEMENTS_BED"
  save_params "$MAF2BED_PARAMS_FILE" "$MAF2BED_PARAMS_STR"
  log "STEP maf2bed END"
else
  log "STEP maf2bed SKIP (output + tham so da khop; dung --force de chay lai)"
fi

# --- Buoc 5: merge-beds (gop voi bulbul_core.bed - dau ra paf2bed cua A2) ---
if [ ! -f "$BULBUL_CORE_BED" ]; then
  echo "run_last.sh: thieu $BULBUL_CORE_BED - chay p02_core_map/run_a2.sh (het buoc paf2bed) truoc buoc merge-beds" >&2
  exit 6
fi
MERGE_PARAMS_FILE="$(params_file_for "$BULBUL_CORE_ALL_BED")"
MERGE_PARAMS_STR="in=[$(input_id "$BULBUL_CORE_BED" "$BULBUL_SHORT_BED")]"
if step_needed "$MERGE_PARAMS_FILE" "$MERGE_PARAMS_STR" "$BULBUL_CORE_ALL_BED" "$MERGE_BEDS_STATS"; then
  log "STEP merge_beds START"
  python -B "$PY" merge-beds --in "$BULBUL_CORE_BED" --in "$BULBUL_SHORT_BED" \
    --out "$BULBUL_CORE_ALL_BED" --stats "$MERGE_BEDS_STATS"
  save_params "$MERGE_PARAMS_FILE" "$MERGE_PARAMS_STR"
  log "STEP merge_beds END"
else
  log "STEP merge_beds SKIP (output + tham so da khop; dung --force de chay lai)"
fi

log "=== run_last.sh xong. Dau ra chinh: $BULBUL_SHORT_BED , $BULBUL_CORE_ALL_BED ==="
