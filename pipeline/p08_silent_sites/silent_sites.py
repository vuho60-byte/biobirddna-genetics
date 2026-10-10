#!/usr/bin/env python3
"""silent_sites.py -- loi X10 / viec B11: phep thu VI TRI CAM cho gene Hox (va
bo gene cong cu phat trien) giua NGUOI va CHAO MAO (Pycnonotus jocosus,
GCA_013400435.1).

Cau hoi (docs/00 Q4, thao luan 2026-10-10): nguoi va chim dung chung cum Hox
va bo gene cong cu. Hai cach giai thich cho hai du doan KHAC NHAU o cac vi tri
KHONG anh huong protein (vi tri codon thu 3 "4-fold degenerate", goi la vi tri
cam):

  H1/H2 (thua ke tu to tien chung ~320 Ma, chon loc giu protein):
         protein giong nhau cao, vi tri cam BAO HOA dot bien -> do giong o vi
         tri cam tien ve muc NGAU NHIEN (~0,25-0,4 tuy thanh phan GC).
  "Ghep/chep gan day" (gene duoc dua tu mot bo gene sang bo gene khac trong
         thoi gian ngan so voi toc do dot bien trung tinh):
         vi tri cam GAN GIONG HET (~1,0) du protein co giong hay khong.

Hai kich ban cho hai con so rat khac nhau; day la mot phep do, khong phai mot
y kien. Nhan toi da cho ket qua: STRONG_INFERENCE (mot cap loai, mot bo gene
chao mao muc scaffold) cho toi khi lap lai voi loai thu hai (ga/se van).

Lenh con (Python 3.11 stdlib; khong numpy/Biopython; khong import module khac
trong pipeline/ -- moi loi tu chua):

  fetch-human  NCBI E-utilities (esearch nuccore + efetch fasta_cds_na)
               -> human_cds.fna + manifest + 1 dong data/registry.csv. CAN MANG.
  compare      human_cds.fna + cds_from_genomic.fna.gz (D01, doc theo luong)
               -> can protein (Needleman-Wunsch/Gotoh, BLOSUM62) -> anh xa ve
               codon -> silent_sites.tsv + alignments.txt (kiem thu cong).
  report       silent_sites.tsv -> report.md: bang theo lop (hox / toolkit /
               control), muc ngau nhien, verdict tung gene, canh bao bat buoc.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from contextlib import contextmanager

# ---------------------------------------------------------------------------
# Tien ich (giong quy uoc pipeline/p04_anomaly, p07_tooth_genes)
# ---------------------------------------------------------------------------


@contextmanager
def smart_open_read(path):
    if path == "-":
        yield sys.stdin
        return
    if str(path).endswith(".gz"):
        f = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        f = open(path, "r", newline="", encoding="utf-8")
    try:
        yield f
    finally:
        f.close()


@contextmanager
def smart_open_write(path):
    if path == "-":
        yield sys.stdout
        return
    f = open(path, "w", newline="", encoding="utf-8")
    try:
        yield f
    finally:
        f.close()


def _fmt_cell(value):
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        return "{:.6g}".format(value)
    return str(value)


def _write_tsv(fileobj, header, rows):
    fileobj.write("\t".join(header) + "\n")
    for row in rows:
        fileobj.write("\t".join(_fmt_cell(row[h]) for h in header) + "\n")


def _read_tsv(path):
    with smart_open_read(path) as f:
        header = None
        rows = []
        for raw in f:
            line = raw.rstrip("\r\n")
            if line == "":
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return header or [], rows


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_fasta(fileobj):
    seq_id = None
    desc = ""
    chunks = []
    for raw in fileobj:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line[0] == ">":
            if seq_id is not None:
                yield seq_id, desc, "".join(chunks)
            parts = line[1:].split(None, 1)
            seq_id = parts[0] if parts else ""
            desc = parts[1] if len(parts) > 1 else ""
            chunks = []
        else:
            if seq_id is None:
                raise SystemExit("FASTA co du lieu truoc dong header '>' dau tien")
            chunks.append(line.strip())
    if seq_id is not None:
        yield seq_id, desc, "".join(chunks)


def write_fasta_record(fileobj, header, seq, width=60):
    fileobj.write(">" + header + "\n")
    for i in range(0, len(seq), width):
        fileobj.write(seq[i : i + width] + "\n")


_HEADER_ATTR_RE = re.compile(r"\[(\w+)=([^\[\]]*)\]")


def header_attrs(desc):
    """'[gene=HOXA1] [protein_id=NP_1.1] ...' -> dict (NCBI cds_from_genomic /
    fasta_cds_na)."""
    return {k: v for k, v in _HEADER_ATTR_RE.findall(desc)}


# ---------------------------------------------------------------------------
# Bang ma, dich, 4-fold degenerate
# ---------------------------------------------------------------------------

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
# Ho codon ma 2 base dau quyet dinh axit amin (vi tri 3 hoan toan thoai hoa):
FOURFOLD_PREFIXES = frozenset(["CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"])


def translate(nt):
    nt = nt.upper()
    return "".join(CODON_TABLE.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))


def is_fourfold_pair(c1, c2):
    """Cap codon (cung cot can) la vi tri 4-fold khi 2 base dau GIONG NHAU va
    thuoc ho 4-fold -> moi thay doi o base 3 deu dong nghia."""
    return c1[:2] == c2[:2] and c1[:2] in FOURFOLD_PREFIXES and set(c1 + c2) <= set("ACGT")


def jukes_cantor(p):
    """Khoang cach JC69 tu p-distance; p >= 0,75 -> inf (bao hoa)."""
    if p <= 0:
        return 0.0
    if p >= 0.75:
        return math.inf
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


# ---------------------------------------------------------------------------
# BLOSUM62 + Needleman-Wunsch/Gotoh (global, affine)
# ---------------------------------------------------------------------------

_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_ROWS = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]
BLOSUM62 = {}
for _i, _a in enumerate(_B62_ORDER):
    for _j, _b in enumerate(_B62_ORDER):
        BLOSUM62[(_a, _b)] = _B62_ROWS[_i][_j]


def b62(a, b):
    """Diem BLOSUM62; ky tu la (X, *, B, Z...) -> -1 (trung tinh nhe)."""
    return BLOSUM62.get((a, b), -1)


def align_global(p1, p2, gap_open=11, gap_extend=1):
    """Can toan cuc 2 protein (Gotoh, gap affine). -> (aln1, aln2, score).
    Chi stdlib; O(n*m) bo nho (3 ma tran int) -- Hox ~250-400 aa: duoi 1 s."""
    n, m = len(p1), len(p2)
    NEG = -10**9
    # M: ket thuc bang khop; X: gap trong p2 (p1 chay); Y: gap trong p1
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        Y[0][j] = -gap_open - (j - 1) * gap_extend
    for i in range(1, n + 1):
        a = p1[i - 1]
        Mi, Xi, Yi = M[i], X[i], Y[i]
        Mp, Xp, Yp = M[i - 1], X[i - 1], Y[i - 1]
        for j in range(1, m + 1):
            s = b62(a, p2[j - 1])
            Mi[j] = max(Mp[j - 1], Xp[j - 1], Yp[j - 1]) + s
            Xi[j] = max(Mp[j] - gap_open, Xp[j] - gap_extend, Yp[j] - gap_open)
            Yi[j] = max(Mi[j - 1] - gap_open, Yi[j - 1] - gap_extend, Xi[j - 1] - gap_open)
    # traceback
    i, j = n, m
    state = max((M[n][m], "M"), (X[n][m], "X"), (Y[n][m], "Y"))[1]
    score = max(M[n][m], X[n][m], Y[n][m])
    a1, a2 = [], []
    while i > 0 or j > 0:
        if state == "M":
            a1.append(p1[i - 1])
            a2.append(p2[j - 1])
            s = b62(p1[i - 1], p2[j - 1])
            target = M[i][j] - s
            i -= 1
            j -= 1
            if i == 0 and j == 0:
                break
            if M[i][j] == target:
                state = "M"
            elif X[i][j] == target:
                state = "X"
            else:
                state = "Y"
        elif state == "X":
            a1.append(p1[i - 1])
            a2.append("-")
            cur = X[i][j]
            i -= 1
            if i == 0 and j == 0:
                break
            if j == 0:
                state = "X"
            elif M[i][j] - gap_open == cur:
                state = "M"
            elif X[i][j] - gap_extend == cur:
                state = "X"
            else:
                state = "Y"
        else:
            a1.append("-")
            a2.append(p2[j - 1])
            cur = Y[i][j]
            j -= 1
            if i == 0 and j == 0:
                break
            if i == 0:
                state = "Y"
            elif M[i][j] - gap_open == cur:
                state = "M"
            elif Y[i][j] - gap_extend == cur:
                state = "Y"
            else:
                state = "X"
    return "".join(reversed(a1)), "".join(reversed(a2)), score


# ---------------------------------------------------------------------------
# Thong ke tren can codon
# ---------------------------------------------------------------------------


def codon_columns(aln1, aln2, cds1, cds2):
    """Tu can protein -> danh sach cot (codon1, codon2, aa1, aa2) chi gom cot
    KHONG gap o ca 2 ben. cds da bo codon dung cuoi."""
    cols = []
    i1 = i2 = 0
    for a, b in zip(aln1, aln2):
        c1 = cds1[3 * i1 : 3 * i1 + 3] if a != "-" else None
        c2 = cds2[3 * i2 : 3 * i2 + 3] if b != "-" else None
        if a != "-":
            i1 += 1
        if b != "-":
            i2 += 1
        if c1 is not None and c2 is not None and len(c1) == 3 and len(c2) == 3:
            cols.append((c1, c2, a, b))
    return cols


def core_mask(cols, window=10, min_ident=0.5):
    """Cot 'loi' = nam trong cua so +/-window cot (khong gap) co ti le aa giong
    >= min_ident. Loai vung can mo ho (duoi Hox ngoai homeodomain, vung lap
    poly-Ala/Pro) khoi thong ke chinh."""
    n = len(cols)
    ident = [1 if c[2] == c[3] else 0 for c in cols]
    prefix = [0]
    for v in ident:
        prefix.append(prefix[-1] + v)
    mask = []
    for i in range(n):
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        frac = (prefix[hi] - prefix[lo]) / (hi - lo)
        mask.append(frac >= min_ident)
    return mask


def site_stats(cols, mask=None):
    """-> dict: n_codon, aa_ident, ident_p1, ident_p2, ident_p3, n_4fold,
    ident_4fold, chance_4fold, jc_4fold, gc3_1, gc3_2."""
    sel = [c for c, m in zip(cols, mask)] if mask is not None else list(cols)
    sel = [c for c, m in zip(cols, mask) if m] if mask is not None else sel
    n = len(sel)
    out = {"n_codon": n, "aa_ident": 0.0, "ident_p1": 0.0, "ident_p2": 0.0, "ident_p3": 0.0,
           "n_4fold": 0, "ident_4fold": 0.0, "chance_4fold": 0.0, "jc_4fold": 0.0, "gc3_1": 0.0, "gc3_2": 0.0}
    if n == 0:
        return out
    aa = p1 = p2 = p3 = 0
    gc1 = gc2 = 0
    ff = 0
    ff_ident = 0
    f1 = {b: 0 for b in "ACGT"}
    f2 = {b: 0 for b in "ACGT"}
    for c1, c2, a, b in sel:
        aa += a == b
        p1 += c1[0] == c2[0]
        p2 += c1[1] == c2[1]
        p3 += c1[2] == c2[2]
        gc1 += c1[2] in "GC"
        gc2 += c2[2] in "GC"
        if is_fourfold_pair(c1, c2):
            ff += 1
            ff_ident += c1[2] == c2[2]
            f1[c1[2]] += 1
            f2[c2[2]] += 1
    out.update({
        "aa_ident": aa / n, "ident_p1": p1 / n, "ident_p2": p2 / n, "ident_p3": p3 / n,
        "gc3_1": gc1 / n, "gc3_2": gc2 / n, "n_4fold": ff,
    })
    if ff:
        ident = ff_ident / ff
        chance = sum((f1[b] / ff) * (f2[b] / ff) for b in "ACGT")
        out.update({"ident_4fold": ident, "chance_4fold": chance, "jc_4fold": jukes_cantor(1.0 - ident)})
    return out


def strip_stop(cds):
    cds = cds.upper()
    if len(cds) >= 3 and CODON_TABLE.get(cds[-3:]) == "*":
        return cds[:-3]
    return cds


def check_cds(cds):
    """-> (ok, ly_do). Yeu cau: chia het 3, >= 30 codon, khong codon dung noi."""
    if len(cds) % 3 != 0:
        return False, "len_not_multiple_of_3"
    if len(cds) < 90:
        return False, "too_short"
    prot = translate(cds)
    if "*" in prot:
        return False, "internal_stop"
    return True, ""


# ---------------------------------------------------------------------------
# Lop gene
# ---------------------------------------------------------------------------

DEFAULT_HOX = [c + str(n) for c in "ABCD" for n in range(1, 14)]  # 52 ten; nguoi co 39
DEFAULT_TOOLKIT = ["PAX6", "SOX9", "RUNX2", "SHH", "BMP4", "TBX5", "TBX4", "PITX1", "FGF8", "GLI3", "HAND2"]
DEFAULT_CONTROL = ["ACTB", "GAPDH", "EEF1A1", "HSPA8", "PGK1", "TPI1", "ENO1", "LDHA", "ALDOA", "RPLP0"]
_HOX_RE = re.compile(r"^HOX[ABCD]\d{1,2}$", re.I)


def gene_class(gene, classes=None):
    g = gene.upper()
    if classes and g in classes:
        return classes[g]
    if _HOX_RE.match(g):
        return "hox"
    if g in DEFAULT_TOOLKIT:
        return "toolkit"
    return "control"


# ---------------------------------------------------------------------------
# fetch-human
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_URLOPEN = urllib.request.urlopen
USER_AGENT = "birdbiodna-silent_sites/0.1 (python-urllib)"


REQUEST_INTERVAL = 0.34  # giay; NCBI gioi han 3 yeu cau/giay khi khong co API key (test dat 0)
_LAST_REQUEST = [0.0]


def _http_get_text(url, timeout=60):
    # Gian cach toi thieu giua 2 yeu cau de khong vuot gioi han NCBI (loi 429).
    wait = REQUEST_INTERVAL - (time.monotonic() - _LAST_REQUEST[0])
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with _URLOPEN(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8")
    finally:
        _LAST_REQUEST[0] = time.monotonic()


def _eutils_url(tool, params, email, api_key):
    q = dict(params)
    q["tool"] = "birdbiodna_silent_sites"
    q["email"] = email
    if api_key:
        q["api_key"] = api_key
    return EUTILS + tool + "?" + urllib.parse.urlencode(q)


def esearch_mrna_ids(gene, organism, email, api_key=None, retmax=20):
    term = '{g}[Gene Name] AND "{o}"[Organism] AND refseq[Filter] AND biomol_mrna[PROP]'.format(g=gene, o=organism)
    url = _eutils_url("esearch.fcgi", {"db": "nuccore", "term": term, "retmode": "json", "retmax": str(retmax)}, email, api_key)
    data = json.loads(_http_get_text(url))
    return list(data.get("esearchresult", {}).get("idlist", [])), term


def efetch_cds_na(ids, email, api_key=None):
    if not ids:
        return ""
    url = _eutils_url("efetch.fcgi", {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta_cds_na", "retmode": "text"}, email, api_key)
    return _http_get_text(url)


MANIFEST_HEADER = ["qseqid", "gene", "class", "organism", "mrna_acc", "protein_id", "length_nt", "n_candidates", "status", "title", "search_term"]


def _acc_from_lcl(seq_id):
    """'lcl|NM_005522.5_cds_NP_005513.1_1' -> 'NM_005522.5'."""
    s = seq_id.split("|", 1)[1] if seq_id.startswith("lcl|") else seq_id
    return s.split("_cds_")[0]


def fetch_human(genes, organism, email, api_key=None, retmax=20, classes=None, log=None):
    fasta_out = []
    manifest = []
    for gene in genes:
        cls = gene_class(gene, classes)
        ids, term = esearch_mrna_ids(gene, organism, email, api_key, retmax)
        if log:
            log("fetch-human: {g} ({c}) esearch -> {n} id".format(g=gene, c=cls, n=len(ids)))
        if not ids:
            manifest.append({"qseqid": "", "gene": gene, "class": cls, "organism": organism, "mrna_acc": "", "protein_id": "",
                             "length_nt": 0, "n_candidates": 0, "status": "not_found", "title": "", "search_term": term})
            continue
        text = efetch_cds_na(ids, email, api_key)
        cands = []
        for sid, desc, seq in iter_fasta(text.splitlines(True)):
            attrs = header_attrs(desc)
            if attrs.get("gene", "").upper() != gene.upper():
                continue
            seq = seq.upper()
            ok, _why = check_cds(strip_stop(seq))
            if not ok:
                continue
            cands.append((_acc_from_lcl(sid), attrs.get("protein_id", ""), attrs.get("protein", ""), seq))
        if not cands:
            manifest.append({"qseqid": "", "gene": gene, "class": cls, "organism": organism, "mrna_acc": "", "protein_id": "",
                             "length_nt": 0, "n_candidates": 0, "status": "no_valid_cds", "title": "", "search_term": term})
            continue
        acc, pid, title, seq = sorted(cands, key=lambda c: (-len(c[3]), c[0]))[0]
        qseqid = "{g}|human|{a}|{p}".format(g=gene.upper(), a=acc, p=pid)
        fasta_out.append((qseqid + " [gene={g}] [protein_id={p}] {t}".format(g=gene.upper(), p=pid, t=title), seq))
        manifest.append({"qseqid": qseqid, "gene": gene.upper(), "class": cls, "organism": organism, "mrna_acc": acc, "protein_id": pid,
                         "length_nt": len(seq), "n_candidates": len(cands), "status": "ok", "title": title, "search_term": term})
    return fasta_out, manifest


def append_registry_row(registry_path, dataset_id, name, source_url, version, local_path, size_bytes, sha, notes, today):
    with open(registry_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row and row[0] == dataset_id:
                return False
    values = {"dataset_id": dataset_id, "name": name, "source_url": source_url, "version_or_release": version,
              "downloaded_on": today, "local_path": local_path, "size_bytes": str(size_bytes), "sha256": sha,
              "license": "NCBI public", "notes": notes}
    with open(registry_path, "a", newline="", encoding="utf-8") as f:
        csv.writer(f, lineterminator="\n").writerow([values.get(h, "") for h in header])
    return True


def load_classes(path):
    if not path:
        return None
    _, rows = _read_tsv(path)
    return {r["gene"].upper(): r["class"] for r in rows if "gene" in r and "class" in r}


def cmd_fetch_human(args):
    if not args.email or "@" not in args.email:
        raise SystemExit("fetch-human: can --email hop le (NCBI E-utilities yeu cau)")
    genes = [g for g in args.genes.split(",") if g]
    classes = load_classes(args.classes)

    def log(msg):
        sys.stderr.write(msg + "\n")

    fasta_out, manifest = fetch_human(genes, args.organism, args.email, args.api_key, args.retmax, classes, log)
    with smart_open_write(args.out) as f:
        for h, s in fasta_out:
            write_fasta_record(f, h, s)
    with smart_open_write(args.manifest) as f:
        _write_tsv(f, MANIFEST_HEADER, manifest)
    n_ok = sum(1 for m in manifest if m["status"] == "ok")
    log("fetch-human: {ok}/{n} gene co CDS; ghi {o} + {m}".format(ok=n_ok, n=len(genes), o=args.out, m=args.manifest))
    if args.registry and args.out != "-" and n_ok:
        sha = sha256_file(args.out)
        today = args.date or _dt.date.today().isoformat()
        accs = ";".join(m["mrna_acc"] for m in manifest if m["status"] == "ok")
        added = append_registry_row(args.registry, "D10-b11-human-cds-" + sha[:8],
                                    "B11 CDS nguoi (RefSeq mRNA) cho phep thu vi tri cam ({o})".format(o=args.organism),
                                    EUTILS + "efetch.fcgi?db=nuccore&rettype=fasta_cds_na", "NCBI RefSeq " + today,
                                    args.out, os.path.getsize(args.out), sha, "accession: " + accs + "; manifest: " + args.manifest, today)
        log("fetch-human: registry " + ("them 1 dong" if added else "da co dong cung sha, khong them"))
    return 0


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

STATS_HEADER = [
    "gene", "class", "human_id", "bulbul_id", "human_len_nt", "bulbul_len_nt", "aln_len", "aln_score",
    "n_codon_all", "aa_ident_all", "ident_p3_all", "n_4fold_all", "ident_4fold_all",
    "n_codon_core", "aa_ident_core", "ident_p1_core", "ident_p2_core", "ident_p3_core",
    "n_4fold_core", "ident_4fold_core", "chance_4fold_core", "jc_4fold_core", "gc3_human", "gc3_bulbul", "status",
]


def _gene_of(seq_id, desc):
    """Ten gene tu header: 'GENE|human|...' hoac '[gene=...]'."""
    attrs = header_attrs(desc)
    if "gene" in attrs:
        return attrs["gene"]
    parts = seq_id.split("|")
    if len(parts) >= 2 and parts[1] == "human":
        return parts[0]
    return seq_id


def load_cds_by_gene(path, wanted=None):
    """FASTA(.gz) -> {GENE_UPPER: (seq_id, cds_upper)} giu CDS DAI NHAT moi
    gene (nhieu isoform). Doc theo luong; wanted (set upper) de bo qua gene
    khong can -> bo nho nho du file CDS chao mao ~40 MB giai nen."""
    best = {}
    with smart_open_read(path) as f:
        for sid, desc, seq in iter_fasta(f):
            g = _gene_of(sid, desc).upper()
            if not g or (wanted is not None and g not in wanted):
                continue
            seq = strip_stop(seq.upper())
            cur = best.get(g)
            if cur is None or len(seq) > len(cur[1]) or (len(seq) == len(cur[1]) and sid < cur[0]):
                best[g] = (sid, seq)
    return best


def compare_pair(gene, cls, hid, hcds, bid, bcds, window=10, min_ident=0.5):
    row = {h: "" for h in STATS_HEADER}
    row.update({"gene": gene, "class": cls, "human_id": hid, "bulbul_id": bid, "human_len_nt": len(hcds), "bulbul_len_nt": len(bcds)})
    for k in STATS_HEADER:
        if k.startswith(("n_", "aa_", "ident_", "chance_", "jc_", "gc3_", "aln_")):
            row[k] = 0
    for who, cds in (("human", hcds), ("bulbul", bcds)):
        ok, why = check_cds(cds)
        if not ok:
            row["status"] = who + "_" + why
            return row, None
    hp, bp = translate(hcds), translate(bcds)
    a1, a2, score = align_global(hp, bp)
    cols = codon_columns(a1, a2, hcds, bcds)
    mask = core_mask(cols, window, min_ident)
    s_all = site_stats(cols)
    s_core = site_stats(cols, mask)
    row.update({
        "aln_len": len(a1), "aln_score": score,
        "n_codon_all": s_all["n_codon"], "aa_ident_all": s_all["aa_ident"], "ident_p3_all": s_all["ident_p3"],
        "n_4fold_all": s_all["n_4fold"], "ident_4fold_all": s_all["ident_4fold"],
        "n_codon_core": s_core["n_codon"], "aa_ident_core": s_core["aa_ident"], "ident_p1_core": s_core["ident_p1"],
        "ident_p2_core": s_core["ident_p2"], "ident_p3_core": s_core["ident_p3"], "n_4fold_core": s_core["n_4fold"],
        "ident_4fold_core": s_core["ident_4fold"], "chance_4fold_core": s_core["chance_4fold"], "jc_4fold_core": s_core["jc_4fold"],
        "gc3_human": s_core["gc3_1"], "gc3_bulbul": s_core["gc3_2"],
        "status": "ok" if s_core["n_4fold"] >= 20 else "few_4fold_sites",
    })
    return row, (a1, a2, mask, cols)


def _alignment_text(gene, hid, bid, a1, a2, mask_cols, width=60):
    """Van ban can protein + dong danh dau loi ('#' = cot loi, '.' = ngoai loi,
    ' ' = gap) de kiem thu cong."""
    a1_m, a2_m, mask, cols = mask_cols
    # mask theo cot khong gap -> dong danh dau theo cot can
    marks = []
    k = 0
    for x, y in zip(a1, a2):
        if x != "-" and y != "-":
            marks.append("#" if mask[k] else ".")
            k += 1
        else:
            marks.append(" ")
    marks = "".join(marks)
    lines = ["## {g}  human={h}  bulbul={b}  aln_len={n}".format(g=gene, h=hid, b=bid, n=len(a1))]
    for i in range(0, len(a1), width):
        lines.append("human  " + a1[i : i + width])
        lines.append("bulbul " + a2[i : i + width])
        lines.append("core   " + marks[i : i + width])
        lines.append("")
    return "\n".join(lines) + "\n"


def cmd_compare(args):
    classes = load_classes(args.classes)
    human = load_cds_by_gene(args.human)
    if not human:
        raise SystemExit("compare: --human khong co CDS nao")
    wanted = set(human)
    bulbul = load_cds_by_gene(args.bulbul, wanted)
    rows = []
    aln_chunks = []
    for gene in sorted(human):
        hid, hcds = human[gene]
        cls = gene_class(gene, classes)
        if gene not in bulbul:
            row = {h: 0 for h in STATS_HEADER}
            row.update({"gene": gene, "class": cls, "human_id": hid, "bulbul_id": "", "human_len_nt": len(hcds), "bulbul_len_nt": 0, "status": "no_bulbul_cds"})
            rows.append(row)
            continue
        bid, bcds = bulbul[gene]
        row, aln = compare_pair(gene, cls, hid, hcds, bid, bcds, args.window, args.min_ident)
        rows.append(row)
        if aln is not None:
            a1, a2, mask, cols = aln
            aln_chunks.append(_alignment_text(gene, hid, bid, a1, a2, (a1, a2, mask, cols)))
    rows.sort(key=lambda r: ({"hox": 0, "toolkit": 1}.get(r["class"], 2), r["gene"]))
    with smart_open_write(args.out) as f:
        _write_tsv(f, STATS_HEADER, rows)
    if args.alignments:
        with smart_open_write(args.alignments) as f:
            f.write("".join(aln_chunks))
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    sys.stderr.write("compare: {n} gene nguoi; {b} co CDS chao mao; {ok} cap du vi tri 4-fold\n".format(
        n=len(rows), b=sum(1 for r in rows if r["bulbul_id"]), ok=n_ok))
    return 0


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

MANDATORY_WARNING = """## Cách đọc — BẮT BUỘC

Phép thử này đo **một** đại lượng: độ giống ở vị trí câm (base thứ 3 của codon 4-fold, nơi mọi thay đổi đều
không đổi axit amin) giữa người và chào mào, so với **mức ngẫu nhiên** tính từ thành phần base của chính hai
trình tự. Hai kịch bản cho hai kỳ vọng được viết ra **trước khi** nhìn số:

| Kịch bản | Protein | Vị trí câm (4-fold) | JC |
|---|---|---|---|
| Thừa kế từ tổ tiên chung ~320 Ma, chọn lọc giữ protein (H1/H2) | giống cao ở vùng lõi | tiến về mức ngẫu nhiên | lớn hoặc bão hòa |
| Chép/ghép gần đây giữa hai bộ gene | bất kỳ | gần 1,0 | gần 0 |

Giới hạn:
1. **Một cặp loài, một bộ gene mức scaffold.** CDS chào mào là mô hình gene dự đoán; exon sai sẽ hạ `aa_ident`
   và kéo `ident_4fold` về ngẫu nhiên một cách giả tạo. Vì vậy thống kê chính tính trên **cột lõi** (cửa sổ
   21 aa có ≥ 50 % aa giống) và chỉ gene có ≥ 20 vị trí 4-fold lõi mới được xét. Lặp lại với gà/sẻ vằn
   trước khi ghi nhãn cao hơn STRONG_INFERENCE.
2. **Thành phần GC khác nhau** giữa hai loài làm mức ngẫu nhiên lệch khỏi 0,25; cột `chance_4fold_core` đã
   tính từ chính dữ liệu, so với cột đó, không so với 0,25.
3. **Chọn lọc lên codon** (codon usage bias, vùng điều hòa trong CDS) có thể giữ vị trí câm trên mức ngẫu nhiên
   một chút; điều đó **không** tạo ra giá trị gần 1,0 trên cả gene.
4. Gene đối chứng (housekeeping) phải cho cùng bức tranh; nếu Hox khác hẳn đối chứng theo hướng "giống hệt
   vị trí câm", đó mới là tín hiệu đáng xem.
5. Kết quả **không** phân biệt được "thừa kế" với một kịch bản thiết kế cố tình làm ngẫu nhiên vị trí câm; nhưng
   kịch bản đó không còn là "ghép/chép" và phải tự nêu dự đoán riêng (đúng tinh thần ô H3 trong `docs/00`).

Cấm viết "đã chứng minh". Không nâng nhãn.
"""


def verdict(row, hi=0.9, band=0.15):
    if row["status"] != "ok":
        return "không xét ({s})".format(s=row["status"])
    ident = float(row["ident_4fold_core"])
    chance = float(row["chance_4fold_core"])
    if ident >= hi:
        return "vị trí câm gần giống hệt → dấu hiệu chép/ghép gần đây, KIỂM LẠI thủ công"
    if abs(ident - chance) <= band:
        return "vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc"
    if ident > chance:
        return "trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây"
    return "dưới ngẫu nhiên → kiểm căn trình tự / thành phần"


def cmd_report(args):
    _, rows = _read_tsv(args.tsv)
    today = args.date or _dt.date.today().isoformat()
    out = ["# B11 — phép thử vị trí câm: Hox và bộ gene công cụ, người so với chào mào", "",
           "Ngày tạo: {d} · công cụ: `pipeline/p08_silent_sites/silent_sites.py report` · đầu vào: `{t}`".format(d=today, t=args.tsv), "",
           "Chữ ký mọi kết quả: **STRONG_INFERENCE tối đa** — đọc mục cuối trước khi trích dẫn.", ""]
    for cls, title in (("hox", "1. Gene Hox"), ("toolkit", "2. Bộ gene công cụ phát triển"), ("control", "3. Gene đối chứng (housekeeping)")):
        sub = [r for r in rows if r["class"] == cls]
        out.append("## " + title)
        out.append("")
        if not sub:
            out.append("_Không có gene nào trong lớp này._")
            out.append("")
            continue
        out.append("| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |")
        out.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in sub:
            if r["status"] in ("no_bulbul_cds",):
                out.append("| {g} | {h} | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |".format(g=r["gene"], h=r["human_id"]))
                continue
            out.append("| {g} | {h} | {b} | {n} | {aa} | {f} | {i} | {c} | {jc} | {v} |".format(
                g=r["gene"], h=r["human_id"], b=r["bulbul_id"], n=r["n_codon_core"],
                aa=_pct(r["aa_ident_core"]), f=r["n_4fold_core"], i=_pct(r["ident_4fold_core"]), c=_pct(r["chance_4fold_core"]),
                jc=r["jc_4fold_core"], v=verdict(r, args.hi, args.band)))
        ok = [r for r in sub if r["status"] == "ok"]
        if ok:
            med_i = _median([float(r["ident_4fold_core"]) for r in ok])
            med_c = _median([float(r["chance_4fold_core"]) for r in ok])
            med_aa = _median([float(r["aa_ident_core"]) for r in ok])
            out.append("")
            out.append("Trung vị lớp này: aa giống {a}, vị trí 4-fold giống {i}, ngẫu nhiên {c} ({n} gene đủ dữ liệu).".format(
                a=_pct(med_aa), i=_pct(med_i), c=_pct(med_c), n=len(ok)))
        out.append("")
    hox_ok = [float(r["ident_4fold_core"]) for r in rows if r["class"] == "hox" and r["status"] == "ok"]
    ctl_ok = [float(r["ident_4fold_core"]) for r in rows if r["class"] == "control" and r["status"] == "ok"]
    out.append("## 4. So Hox với đối chứng")
    out.append("")
    if hox_ok and ctl_ok:
        out.append("Trung vị giống 4-fold: Hox {h} so với đối chứng {c}. Nếu hai số này cùng nằm gần mức ngẫu nhiên, Hox không khác gì phần còn lại của bộ gene ở vị trí câm.".format(
            h=_pct(_median(hox_ok)), c=_pct(_median(ctl_ok))))
    else:
        out.append("_Chưa đủ gene ở cả hai lớp để so._")
    out.append("")
    out.append(MANDATORY_WARNING)
    with smart_open_write(args.out) as f:
        f.write("\n".join(out))
    return 0


def _pct(x):
    try:
        return "{:.1%}".format(float(x))
    except (TypeError, ValueError):
        return str(x)


def _median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return 0.0
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser():
    p = argparse.ArgumentParser(prog="silent_sites.py", description="B11: phép thử vị trí câm Hox/công cụ, người so với chào mào (xem README.md).")
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch-human", help="NCBI esearch nuccore + efetch fasta_cds_na -> human_cds.fna (+manifest, +registry). CẦN MẠNG.")
    f.add_argument("--genes", default=",".join(["HOX" + x for x in DEFAULT_HOX] + DEFAULT_TOOLKIT + DEFAULT_CONTROL),
                   help="mặc định: 52 tên HOX (người có 39, tên thiếu sẽ not_found), 11 gene công cụ, 10 gene đối chứng")
    f.add_argument("--organism", default="Homo sapiens")
    f.add_argument("--email", required=True)
    f.add_argument("--api-key", default=None)
    f.add_argument("--retmax", type=int, default=20)
    f.add_argument("--classes", default=None, help="TSV gene<TAB>class ghi đè phân lớp")
    f.add_argument("--out", required=True)
    f.add_argument("--manifest", required=True)
    f.add_argument("--registry", default=None)
    f.add_argument("--date", default=None)
    f.set_defaults(func=cmd_fetch_human)

    c = sub.add_parser("compare", help="human_cds.fna + CDS chào mào (.fna.gz) -> silent_sites.tsv (+alignments.txt)")
    c.add_argument("--human", required=True)
    c.add_argument("--bulbul", required=True, help="data/ncbi/GCA_013400435.1/..._cds_from_genomic.fna.gz")
    c.add_argument("--classes", default=None)
    c.add_argument("--window", type=int, default=10, help="nửa cửa sổ (aa) để xác định cột lõi")
    c.add_argument("--min-ident", type=float, default=0.5, help="tỉ lệ aa giống tối thiểu trong cửa sổ để là cột lõi")
    c.add_argument("--out", required=True)
    c.add_argument("--alignments", default=None, help="file văn bản căn protein + dòng đánh dấu lõi")
    c.set_defaults(func=cmd_compare)

    r = sub.add_parser("report", help="silent_sites.tsv -> report.md")
    r.add_argument("--tsv", required=True)
    r.add_argument("--hi", type=float, default=0.9, help="ngưỡng 'gần giống hệt' cho vị trí 4-fold")
    r.add_argument("--band", type=float, default=0.15, help="độ rộng quanh mức ngẫu nhiên coi là 'ở mức ngẫu nhiên'")
    r.add_argument("--date", default=None)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_report)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
