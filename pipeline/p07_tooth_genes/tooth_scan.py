#!/usr/bin/env python3
"""tooth_scan.py -- loi X9 / viec B5 (docs/05): quet MANH VO GENE RANG (va GULO)
trong bo gene chao mao (Pycnonotus jocosus, GCA_013400435.1).

Cau hoi khoa hoc (docs/00 Q4, tong hop research/synthesis/11): neu chim hien
dai chia se MOT to tien chung da mat rang khoang hoa (~116 Ma, Meredith 2014
-- R12 dong #23, CONFIRMED), thi 6 gene men/nga rang ENAM, AMBN, AMELX, AMTN,
DSPP, MMP20 phai con lai trong bo gene chao mao duoi dang GIA GENE (pseudogene)
mang ton thuong bat hoat (codon dung som, lech khung doc) o vi tri KHOP voi
cac bo chim khac. Gia thuyet H3 (bo gene duoc thiet ke rieng) khong du doan
mot bo gene "sach" phai mang gene hong vo dung cung kieu loi voi ho hang.

Phuong phap chuan (tong hop 07 muc 5): TBLASTN voi protein ortholog cua
NGOAI NHOM (ca sau / rua -- D06 docs/02) lam truy van, can len bo gene chao
mao; roi KIEM THU CONG dich khung va codon dung. Chua ai lam o bo Se.

Lenh con (Python 3.11 stdlib, khong numpy/Biopython, khong import module khac
trong pipeline/ -- moi loi tu chua, dung quy uoc cua repo):

  fetch-queries  NCBI E-utilities (esearch + efetch) -> queries.faa + manifest
                 + 1 dong data/registry.csv (sha256). CAN MANG. Khong chay o
                 buoc viet code nay.
  parse          bang tblastn (-outfmt 6 voi cot co dinh OUTFMT) -> hits.tsv,
                 loci.tsv (chuoi HSP thanh locus), lesions.tsv (codon dung,
                 lech khung ung vien, ranh exon kieu intron).
  extract        loci.tsv + genome .fna.gz (doc theo luong) -> loci.fa (+pad)
                 va ban dich 3 khung doc de KIEM THU CONG.
  translate      FASTA (nucleotide) -> ban dich 3 khung, danh dau '*'.
  report         loci.tsv + lesions.tsv (+manifest) -> report.md co PHAN
                 CANH BAO BAT BUOC va thu tu kha nang giai thich.

QUAN TRONG: cong cu nay SANG LOC TIN HIEU, khong KET LUAN. Moi dong ket qua
toi da la STRONG_INFERENCE cho toi khi (a) doi chung duong dat, (b) kiem thu
cong vung trich, va (c) doi chieu vi tri ton thuong voi bang phu luc Meredith
2014 (chua co trong repo). "Khong tim thay" KHONG phai "khong ton tai" (lap
rap muc scaffold, N50 ~218 kb; gia gene 116 Ma phan ky manh). Doc README.md.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from contextlib import contextmanager

# ---------------------------------------------------------------------------
# Tien ich dung chung (giong quy uoc pipeline/p04_anomaly/composition_scan.py)
# ---------------------------------------------------------------------------


@contextmanager
def smart_open_read(path):
    """Mo doc text: '-' = stdin; '.gz' = gzip text. newline='' de tu xu ly
    CRLF/LF o noi goi (rstrip('\\r\\n'))."""
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
    """Mo ghi text: '-' = stdout. newline='' de ghi LF thuan tuy."""
    if path == "-":
        yield sys.stdout
        return
    f = open(path, "w", newline="", encoding="utf-8")
    try:
        yield f
    finally:
        f.close()


def _fmt_cell(value):
    """float -> 6 chu so co nghia; con lai -> str() (so nguyen giu nguyen)."""
    if isinstance(value, float):
        return "{:.6g}".format(value)
    return str(value)


def _write_tsv(fileobj, header, rows):
    fileobj.write("\t".join(header) + "\n")
    for row in rows:
        fileobj.write("\t".join(_fmt_cell(row[h]) for h in header) + "\n")


def _read_tsv(path):
    """Doc TSV co header -> (header, [dict,...]). Dong rong bi bo."""
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


# ---------------------------------------------------------------------------
# FASTA + dich ma
# ---------------------------------------------------------------------------


def iter_fasta(fileobj):
    """Sinh (seq_id, mo_ta, trinh_tu) tung record theo luong; chi giu 1 record
    trong bo nho. seq_id = token dau tien sau '>'."""
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


_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMP = str.maketrans("ACGTacgtNnRYKMSWBDHVrykmswbdhv", "TGCAtgcaNnYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def translate(nt, frame=0):
    """Dich tu vi tri `frame` (0/1/2), bang ma chuan; codon co ky tu ngoai
    A/C/G/T -> 'X'; codon dung -> '*'. Khong dung o codon dung."""
    nt = nt.upper()
    out = []
    for i in range(frame, len(nt) - 2, 3):
        out.append(CODON_TABLE.get(nt[i : i + 3], "X"))
    return "".join(out)


# ---------------------------------------------------------------------------
# Dinh danh truy van: 'GENE|ROLE|ACCESSION' (fetch-queries tao) hoac tu do
# ---------------------------------------------------------------------------

ROLE_TARGET = "target"
ROLE_CONTROL = "control"

DEFAULT_TARGET_GENES = ["ENAM", "AMBN", "AMELX", "AMTN", "DSPP", "MMP20", "GULO"]
DEFAULT_CONTROL_GENES = ["MEPE", "SPARC"]


def split_qseqid(qseqid, manifest=None):
    """-> (gene, role, accession). Uu tien manifest (qseqid -> gene, role);
    sau do dang 'GENE|ROLE|ACC'; cuoi cung gene = qseqid, role = target."""
    if manifest and qseqid in manifest:
        m = manifest[qseqid]
        return m.get("gene", qseqid), m.get("role", ROLE_TARGET), m.get("accession", "")
    parts = qseqid.split("|")
    if len(parts) >= 3 and parts[1] in (ROLE_TARGET, ROLE_CONTROL):
        return parts[0], parts[1], "|".join(parts[2:])
    return qseqid, ROLE_TARGET, ""


def load_manifest(path):
    """queries_manifest.tsv -> {qseqid: row}. Tra {} neu path None."""
    if not path:
        return {}
    _, rows = _read_tsv(path)
    return {r["qseqid"]: r for r in rows if "qseqid" in r}


# ---------------------------------------------------------------------------
# fetch-queries: NCBI E-utilities (esearch -> efetch FASTA)
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_URLOPEN = urllib.request.urlopen  # test monkeypatch diem nay
USER_AGENT = "birdbiodna-tooth_scan/0.1 (python-urllib)"


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
    q["tool"] = "birdbiodna_tooth_scan"
    q["email"] = email
    if api_key:
        q["api_key"] = api_key
    return EUTILS + tool + "?" + urllib.parse.urlencode(q)


def esearch_protein_ids(gene, organism, email, api_key=None, refseq_only=True, retmax=20):
    term = '{g}[Gene Name] AND "{o}"[Organism]'.format(g=gene, o=organism)
    if refseq_only:
        term += " AND refseq[Filter]"
    url = _eutils_url(
        "esearch.fcgi",
        {"db": "protein", "term": term, "retmode": "json", "retmax": str(retmax)},
        email,
        api_key,
    )
    data = json.loads(_http_get_text(url))
    return list(data.get("esearchresult", {}).get("idlist", [])), term


def efetch_protein_fasta(ids, email, api_key=None):
    if not ids:
        return ""
    url = _eutils_url(
        "efetch.fcgi",
        {"db": "protein", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
        email,
        api_key,
    )
    return _http_get_text(url)


def _pick_longest(records):
    """records: [(acc, desc, seq)] -> record dai nhat (hoa: accession nho hon
    theo thu tu chuoi, de ket qua on dinh)."""
    return sorted(records, key=lambda r: (-len(r[2]), r[0]))[0]


def fetch_queries(genes, roles, organism, email, api_key=None, refseq_only=True, retmax=20, log=None):
    """-> (fasta_text, manifest_rows). Moi gene: esearch -> efetch -> giu
    isoform DAI NHAT. Gene khong co ban ghi -> dong manifest status=not_found
    (khong crash: thieu ortholog o ngoai nhom la THONG TIN, khong phai loi)."""
    fasta_out = []
    manifest = []
    for gene in genes:
        role = roles.get(gene, ROLE_TARGET)
        ids, term = esearch_protein_ids(gene, organism, email, api_key, refseq_only, retmax)
        if log:
            log("fetch-queries: {g} ({r}) esearch -> {n} id | {t}".format(g=gene, r=role, n=len(ids), t=term))
        if not ids:
            manifest.append(
                {
                    "qseqid": "",
                    "gene": gene,
                    "role": role,
                    "accession": "",
                    "organism": organism,
                    "length_aa": 0,
                    "n_candidates": 0,
                    "status": "not_found",
                    "title": "",
                    "search_term": term,
                }
            )
            continue
        text = efetch_protein_fasta(ids, email, api_key)
        records = []
        for acc, desc, seq in iter_fasta(text.splitlines(True)):
            if seq:
                records.append((acc, desc, seq.upper()))
        if not records:
            manifest.append(
                {
                    "qseqid": "",
                    "gene": gene,
                    "role": role,
                    "accession": "",
                    "organism": organism,
                    "length_aa": 0,
                    "n_candidates": len(ids),
                    "status": "efetch_empty",
                    "title": "",
                    "search_term": term,
                }
            )
            continue
        acc, desc, seq = _pick_longest(records)
        qseqid = "{g}|{r}|{a}".format(g=gene, r=role, a=acc)
        fasta_out.append((qseqid + " " + organism + " " + desc, seq))
        manifest.append(
            {
                "qseqid": qseqid,
                "gene": gene,
                "role": role,
                "accession": acc,
                "organism": organism,
                "length_aa": len(seq),
                "n_candidates": len(records),
                "status": "ok",
                "title": desc,
                "search_term": term,
            }
        )
    return fasta_out, manifest


MANIFEST_HEADER = [
    "qseqid", "gene", "role", "accession", "organism", "length_aa",
    "n_candidates", "status", "title", "search_term",
]


def append_registry_row(registry_path, dataset_id, name, source_url, version, local_path, size_bytes, sha, notes, today):
    """Them 1 dong vao data/registry.csv theo DUNG header hien co (doc header
    truoc, dien theo ten cot). Khong bao gio xoa/ghi de dong cu. Neu da co
    dong cung dataset_id -> khong them (idempotent) -> tra False."""
    with open(registry_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row and row[0] == dataset_id:
                return False
    values = {
        "dataset_id": dataset_id,
        "name": name,
        "source_url": source_url,
        "version_or_release": version,
        "downloaded_on": today,
        "local_path": local_path,
        "size_bytes": str(size_bytes),
        "sha256": sha,
        "license": "NCBI public",
        "notes": notes,
    }
    with open(registry_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow([values.get(h, "") for h in header])
    return True


def cmd_fetch_queries(args):
    targets = [g for g in args.genes.split(",") if g]
    controls = [g for g in args.controls.split(",") if g] if args.controls else []
    roles = {g: ROLE_TARGET for g in targets}
    roles.update({g: ROLE_CONTROL for g in controls})
    genes = targets + [g for g in controls if g not in targets]
    if not args.email or "@" not in args.email:
        raise SystemExit("fetch-queries: can --email hop le (NCBI E-utilities yeu cau)")

    def log(msg):
        sys.stderr.write(msg + "\n")

    fasta_out, manifest = fetch_queries(
        genes, roles, args.organism, args.email, args.api_key,
        refseq_only=not args.all_proteins, retmax=args.retmax, log=log,
    )
    with smart_open_write(args.out) as f:
        for header, seq in fasta_out:
            write_fasta_record(f, header, seq)
    with smart_open_write(args.manifest) as f:
        _write_tsv(f, MANIFEST_HEADER, manifest)
    n_ok = sum(1 for m in manifest if m["status"] == "ok")
    log("fetch-queries: {ok}/{n} gene co protein; ghi {o} + {m}".format(ok=n_ok, n=len(genes), o=args.out, m=args.manifest))
    if args.registry and args.out != "-" and n_ok:
        sha = sha256_file(args.out)
        size = os.path.getsize(args.out)
        today = args.date or _dt.date.today().isoformat()
        accs = ";".join(m["accession"] for m in manifest if m["status"] == "ok")
        added = append_registry_row(
            args.registry,
            dataset_id="D09-b5-queries-" + sha[:8],
            name="B5 truy van protein ngoai nhom ({o})".format(o=args.organism),
            source_url=EUTILS + "efetch.fcgi?db=protein&rettype=fasta",
            version="NCBI Protein " + today,
            local_path=args.out,
            size_bytes=size,
            sha=sha,
            notes="accession: " + accs + "; manifest: " + args.manifest,
            today=today,
        )
        log("fetch-queries: registry " + ("them 1 dong" if added else "da co dong cung sha, khong them"))
    return 0


# ---------------------------------------------------------------------------
# parse: bang tblastn -> HSP -> locus -> ton thuong
# ---------------------------------------------------------------------------

OUTFMT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sframe", "qseq", "sseq",
]
OUTFMT = "6 " + " ".join(OUTFMT_COLUMNS)
_INT_COLS = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen", "sframe"}
_FLOAT_COLS = {"pident", "evalue", "bitscore"}


def read_hsps(path, max_evalue=None):
    """Doc bang tblastn (-outfmt OUTFMT). Dong '#' va dong rong bi bo. Loi so
    cot -> SystemExit ghi ro dong loi (khong doan cot)."""
    hsps = []
    with smart_open_read(path) as f:
        for n, raw in enumerate(f, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != len(OUTFMT_COLUMNS):
                raise SystemExit(
                    "parse: dong {n} co {k} cot, can {m} (chay tblastn voi -outfmt \"{o}\")".format(
                        n=n, k=len(fields), m=len(OUTFMT_COLUMNS), o=OUTFMT
                    )
                )
            h = dict(zip(OUTFMT_COLUMNS, fields))
            for c in _INT_COLS:
                h[c] = int(h[c])
            for c in _FLOAT_COLS:
                h[c] = float(h[c])
            if max_evalue is not None and h["evalue"] > max_evalue:
                continue
            h["strand"] = "+" if h["sframe"] > 0 else "-"
            h["s_lo"] = min(h["sstart"], h["send"])
            h["s_hi"] = max(h["sstart"], h["send"])
            hsps.append(h)
    return hsps


def chain_loci(hsps, max_locus_gap):
    """Nhom HSP theo (qseqid, sseqid, strand), sap theo toa do subject, noi
    thanh locus khi khoang cach subject giua HSP ke nhau <= max_locus_gap
    (nt). -> [ {qseqid, sseqid, strand, hsps:[...]} ] theo thu tu on dinh."""
    groups = {}
    for h in hsps:
        groups.setdefault((h["qseqid"], h["sseqid"], h["strand"]), []).append(h)
    loci = []
    for key in sorted(groups):
        items = sorted(groups[key], key=lambda h: (h["s_lo"], h["s_hi"]))
        current = [items[0]]
        for h in items[1:]:
            prev_hi = max(x["s_hi"] for x in current)
            if h["s_lo"] - prev_hi <= max_locus_gap:
                current.append(h)
            else:
                loci.append({"qseqid": key[0], "sseqid": key[1], "strand": key[2], "hsps": current})
                current = [h]
        loci.append({"qseqid": key[0], "sseqid": key[1], "strand": key[2], "hsps": current})
    return loci


def _union_len(intervals):
    """Tong do dai hop cac doan [a,b] (dong, so nguyen)."""
    total = 0
    cur_a = cur_b = None
    for a, b in sorted(intervals):
        if cur_a is None:
            cur_a, cur_b = a, b
        elif a <= cur_b + 1:
            cur_b = max(cur_b, b)
        else:
            total += cur_b - cur_a + 1
            cur_a, cur_b = a, b
    if cur_a is not None:
        total += cur_b - cur_a + 1
    return total


def _owner_by_qpos(hsps):
    """Vi tri truy van (aa) -> chi so HSP 'chu so huu' = HSP phu vi tri do co
    bitscore cao nhat (hoa: chi so nho hon). Dung de KHONG dem '*' nam trong
    phan HSP yeu mo rong chong lan sang khung sai (hien tuong that quan sat
    duoc tren du lieu tong hop: HSP thu 2 sau lech khung keo dai nguoc 8 aa
    vao vung HSP 1 va sinh '*' gia)."""
    owner = {}
    for idx, h in enumerate(hsps):
        for q in range(h["qstart"], h["qend"] + 1):
            cur = owner.get(q)
            if cur is None or h["bitscore"] > hsps[cur]["bitscore"]:
                owner[q] = idx
    return owner


def _walk_alignment(h):
    """Duyet cot can cua 1 HSP -> list (col, qpos, s_nt, qch, sch).
    qpos = vi tri aa truy van cua cot (cot chen trong subject -> qpos cua aa
    truy van lien truoc). s_nt = toa do nucleotide (1-based, tren mach +)
    cua base DAU codon tuong ung voi cot do (cot gap trong subject -> None)."""
    qpos = h["qstart"] - 1
    s_aa = 0  # so aa subject da di qua
    out = []
    sign = 1 if h["strand"] == "+" else -1
    for col, (qch, sch) in enumerate(zip(h["qseq"], h["sseq"])):
        if qch != "-":
            qpos += 1
        if sch != "-":
            if h["strand"] == "+":
                s_nt = h["sstart"] + 3 * s_aa
            else:
                s_nt = h["sstart"] - 3 * s_aa - 2  # base thap nhat cua codon
            s_aa += 1
        else:
            s_nt = None
        out.append((col, qpos, s_nt, qch, sch))
    return out


def _flank_supported(cols, col_idx, flank, min_ident):
    """Codon dung tai cot col_idx co 'cho dua' khong: moi ben co >= flank cot
    (ca 2 ky tu khong gap) voi ti le khop (qch == sch) >= min_ident."""
    def side(rng):
        taken = 0
        ident = 0
        for c in rng:
            _, _, _, qch, sch = cols[c]
            if qch == "-" or sch == "-":
                continue
            taken += 1
            if qch == sch:
                ident += 1
            if taken == flank:
                break
        return taken >= flank and (ident / taken) >= min_ident

    left = side(range(col_idx - 1, -1, -1))
    right = side(range(col_idx + 1, len(cols)))
    return left and right


def analyze_locus(locus, flank=10, flank_min_ident=0.5, max_q_gap=20, max_fs_gap=40, min_intron=40):
    """-> (summary dict, [lesion dicts]). Xem README muc 'parse' de biet dinh
    nghia tung cot. Khong ket luan gi: chi dem va liet ke."""
    hsps = locus["hsps"]
    qlen = hsps[0]["qlen"]
    qcov = _union_len([(h["qstart"], h["qend"]) for h in hsps]) / qlen if qlen else 0.0
    owner = _owner_by_qpos(hsps)
    lesions = []
    n_stop_raw = 0
    n_stop_supported = 0
    for idx, h in enumerate(hsps):
        cols = _walk_alignment(h)
        for col, qpos, s_nt, qch, sch in cols:
            if sch != "*":
                continue
            n_stop_raw += 1
            owned = owner.get(qpos) == idx
            supported = owned and _flank_supported(cols, col, flank, flank_min_ident)
            if supported:
                n_stop_supported += 1
            lesions.append(
                {
                    "qseqid": locus["qseqid"],
                    "sseqid": locus["sseqid"],
                    "strand": locus["strand"],
                    "type": "stop_codon",
                    "qpos": qpos,
                    "s_nt": s_nt if s_nt is not None else "",
                    "frame_from": h["sframe"],
                    "frame_to": h["sframe"],
                    "gap_nt": 0,
                    "owned": int(owned),
                    "supported": int(supported),
                    "hsp_bitscore": h["bitscore"],
                    "note": "query aa " + str(qpos) + ("" if owned else " (nam trong HSP yeu chong lan, bo qua)"),
                }
            )
    # Lech khung / ranh exon: xet cap HSP ke nhau theo thu tu TRUY VAN
    by_q = sorted(range(len(hsps)), key=lambda i: (hsps[i]["qstart"], hsps[i]["qend"]))
    n_fs = 0
    n_intron = 0
    n_inframe = 0
    for a, b in zip(by_q, by_q[1:]):
        h1, h2 = hsps[a], hsps[b]
        q_gap = h2["qstart"] - h1["qend"]
        if h1["strand"] == "+":
            s_gap = h2["sstart"] - h1["send"]
        else:
            s_gap = h1["send"] - h2["sstart"]
        same_frame = h1["sframe"] == h2["sframe"]
        # Phan loai theo KHOANG CACH SUBJECT truoc (quan sat tren du lieu tong
        # hop: HSP mo rong vao intron co the chong lan tren truy van toi ~50 aa,
        # nen khong dung chong lan truy van de loai ranh exon):
        if s_gap > min_intron:
            kind = "intron_like"  # doi khung giua 2 exon la binh thuong
            n_intron += 1
        elif same_frame:
            kind = "in_frame_break"  # cung khung, dut ngan: indel boi so 3 / vung kho can
            n_inframe += 1
        elif abs(q_gap) <= max_q_gap and s_gap <= max_fs_gap:
            kind = "frameshift"
            n_fs += 1
        else:
            kind = "unclassified"  # khac khung nhung chong lan/gap bat thuong: de nguoi xem
        if kind in ("frameshift", "intron_like", "in_frame_break", "unclassified"):
            lesions.append(
                {
                    "qseqid": locus["qseqid"],
                    "sseqid": locus["sseqid"],
                    "strand": locus["strand"],
                    "type": kind,
                    "qpos": h1["qend"],
                    "s_nt": h1["send"],
                    "frame_from": h1["sframe"],
                    "frame_to": h2["sframe"],
                    "gap_nt": s_gap,
                    "owned": 1,
                    "supported": int(kind == "frameshift"),
                    "hsp_bitscore": min(h1["bitscore"], h2["bitscore"]),
                    "note": "q{a}->q{b}, subject gap {g} nt, khung {f1}->{f2}".format(
                        a=h1["qend"], b=h2["qstart"], g=s_gap, f1=h1["sframe"], f2=h2["sframe"]
                    ),
                }
            )
    summary = {
        "qseqid": locus["qseqid"],
        "sseqid": locus["sseqid"],
        "strand": locus["strand"],
        "s_start": min(h["s_lo"] for h in hsps),
        "s_end": max(h["s_hi"] for h in hsps),
        "n_hsp": len(hsps),
        "qlen": qlen,
        "qcov": qcov,
        "best_evalue": min(h["evalue"] for h in hsps),
        "max_pident": max(h["pident"] for h in hsps),
        "sum_bitscore": sum(h["bitscore"] for h in hsps),
        "frames": ",".join(str(h["sframe"]) for h in sorted(hsps, key=lambda x: x["qstart"])),
        "n_stop_raw": n_stop_raw,
        "n_stop_supported": n_stop_supported,
        "n_frameshift": n_fs,
        "n_intron_like": n_intron,
        "n_in_frame_break": n_inframe,
    }
    return summary, lesions


LOCI_HEADER = [
    "gene", "role", "accession", "qseqid", "sseqid", "strand", "s_start", "s_end", "n_hsp", "qlen",
    "qcov", "best_evalue", "max_pident", "sum_bitscore", "frames", "n_stop_raw",
    "n_stop_supported", "n_frameshift", "n_intron_like", "n_in_frame_break",
]
LESIONS_HEADER = [
    "gene", "role", "qseqid", "sseqid", "strand", "type", "qpos", "s_nt", "frame_from", "frame_to",
    "gap_nt", "owned", "supported", "hsp_bitscore", "note",
]
HITS_HEADER = ["gene", "role"] + OUTFMT_COLUMNS[:-2]  # bo qseq/sseq cho gon


def cmd_parse(args):
    manifest = load_manifest(args.manifest)
    hsps = read_hsps(args.tblastn, max_evalue=args.max_evalue)
    loci = chain_loci(hsps, args.max_locus_gap)
    loci_rows = []
    lesion_rows = []
    for locus in loci:
        summary, lesions = analyze_locus(
            locus, flank=args.flank, flank_min_ident=args.flank_min_ident,
            max_q_gap=args.max_q_gap, max_fs_gap=args.max_fs_gap, min_intron=args.min_intron,
        )
        gene, role, acc = split_qseqid(locus["qseqid"], manifest)
        summary.update({"gene": gene, "role": role, "accession": acc})
        loci_rows.append(summary)
        for les in lesions:
            les.update({"gene": gene, "role": role})
            lesion_rows.append(les)
    loci_rows.sort(key=lambda r: (r["role"] != ROLE_TARGET, r["gene"], -r["sum_bitscore"]))
    lesion_rows.sort(key=lambda r: (r["gene"], r["sseqid"], r["qpos"]))
    with smart_open_write(args.out_loci) as f:
        _write_tsv(f, LOCI_HEADER, loci_rows)
    with smart_open_write(args.out_lesions) as f:
        _write_tsv(f, LESIONS_HEADER, lesion_rows)
    if args.out_hits:
        rows = []
        for h in hsps:
            gene, role, _ = split_qseqid(h["qseqid"], manifest)
            r = dict(h)
            r.update({"gene": gene, "role": role})
            rows.append(r)
        with smart_open_write(args.out_hits) as f:
            _write_tsv(f, HITS_HEADER, rows)
    sys.stderr.write(
        "parse: {h} HSP (evalue <= {e}) -> {l} locus; {s} codon dung co cho dua, {fs} lech khung ung vien\n".format(
            h=len(hsps), e=args.max_evalue, l=len(loci_rows),
            s=sum(r["n_stop_supported"] for r in loci_rows), fs=sum(r["n_frameshift"] for r in loci_rows),
        )
    )
    return 0


# ---------------------------------------------------------------------------
# extract / translate: trich vung locus de kiem thu cong
# ---------------------------------------------------------------------------


def three_frame_report(name, seq):
    """Van ban 3 khung dich cua 1 trinh tu (da o huong doc gene): moi khung 1
    dong protein + danh sach vi tri '*' (toa do nt 1-based trong doan)."""
    lines = ["## " + name, "chieu dai: {n} nt".format(n=len(seq))]
    for fr in range(3):
        prot = translate(seq, fr)
        stops = [fr + 3 * i + 1 for i, aa in enumerate(prot) if aa == "*"]
        lines.append("khung +{f}: {k} codon dung tai nt {s}".format(
            f=fr + 1, k=len(stops), s=",".join(str(x) for x in stops) if stops else "-"))
        lines.append(prot)
    return "\n".join(lines) + "\n"


def cmd_extract(args):
    _, loci = _read_tsv(args.loci)
    wanted = {}
    for r in loci:
        if args.role and r.get("role") != args.role:
            continue
        wanted.setdefault(r["sseqid"], []).append(r)
    if not wanted:
        raise SystemExit("extract: loci.tsv khong co dong nao (sau loc --role)")
    n_written = 0
    with smart_open_read(args.fna) as fin, smart_open_write(args.out) as fout:
        ftr = open(args.translate, "w", newline="", encoding="utf-8") if args.translate else None
        try:
            for sid, _desc, seq in iter_fasta(fin):
                if sid not in wanted:
                    continue
                for r in wanted[sid]:
                    s0 = max(1, int(r["s_start"]) - args.pad)
                    s1 = min(len(seq), int(r["s_end"]) + args.pad)
                    sub = seq[s0 - 1 : s1].upper()
                    if r["strand"] == "-":
                        sub = revcomp(sub)
                    name = "{g}|{s}:{a}-{b}({st})".format(g=r["gene"], s=sid, a=s0, b=s1, st=r["strand"])
                    write_fasta_record(fout, name + " pad={p} qseqid={q}".format(p=args.pad, q=r["qseqid"]), sub)
                    if ftr:
                        ftr.write(three_frame_report(name, sub))
                    n_written += 1
        finally:
            if ftr:
                ftr.close()
    missing = sum(len(v) for v in wanted.values()) - n_written
    sys.stderr.write("extract: ghi {n} vung; {m} locus khong tim thay scaffold trong --fna\n".format(n=n_written, m=missing))
    return 0 if missing == 0 else 1


def cmd_translate(args):
    with smart_open_read(args.fasta) as fin, smart_open_write(args.out) as fout:
        for sid, _desc, seq in iter_fasta(fin):
            fout.write(three_frame_report(sid, seq.upper()))
    return 0


# ---------------------------------------------------------------------------
# report: Markdown + canh bao bat buoc
# ---------------------------------------------------------------------------

MANDATORY_WARNING = """## Cách đọc kết quả — BẮT BUỘC đọc trước khi dùng bất kỳ dòng nào ở trên

Mọi con số trong báo cáo này là **tín hiệu cần kiểm tiếp**, không phải kết luận. Nhãn tối đa cho
bất kỳ phát biểu nào rút từ đây là **STRONG_INFERENCE** (một phương pháp, một bộ gene), cho tới khi
đủ ba điều kiện: (a) đối chứng dương đạt, (b) vùng trích đã kiểm thủ công bằng `extract`/`translate`,
(c) vị trí tổn thương đã đối chiếu với bảng phụ lục Meredith 2014 *Science* 346:1254390 (chưa có trong repo).

**Thứ tự khả năng giải thích cho một codon dừng / lệch khung** — chỉ chuyển mục sau khi loại mục trước:

1. **Lỗi lắp ráp hoặc lỗi gọi base** — bộ gene chào mào ở mức scaffold (N50 ~218 kb); một thay đổi
   1 nt đơn lẻ có thể là lỗi consensus. Kiểm: có đọc được cùng tổn thương ở loài sẻ họ hàng gần không.
2. **Giả tạo do căn trình tự** — HSP của tblastn mở rộng vào vùng không mã hóa (intron, đầu đoạn) với
   độ giống thấp sinh `*` giả; cột `supported` đã lọc một phần (cánh hai bên ≥ `--flank` aa khớp
   ≥ `--flank-min-ident`), nhưng **không thay được kiểm thủ công**.
3. **Ranh exon thật** (`intron_like`) — đổi khung đọc giữa hai exon là bình thường, **không phải** tổn thương.
4. **Tổn thương bất hoạt thật** — và chỉ khi vị trí **khớp** với các bộ chim khác (Meredith 2014) mới
   được gọi là "cùng kiểu hỏng" (dự đoán của H1/H2). Tổn thương riêng chào mào vẫn tương thích với H1
   (giả gene tích thêm đột biến sau khi hỏng) — không phải bằng chứng chống H1.

**"Không tìm thấy" ≠ "không tồn tại".** Gene răng ở chim là giả gene ~116 Ma, phân kỳ mạnh với cá sấu/rùa;
mảnh có thể nằm ở khoảng trống lắp ráp hoặc dưới ngưỡng e-value. Kết quả âm **không kết luận** được gì về H1/H2/H3.

**Đối chứng dương không đạt** (gene đối chứng có codon dừng/lệch khung hoặc độ phủ thấp) → toàn bộ dòng gene
mục tiêu **không diễn giải được** trong lần chạy này: phải sửa tham số hoặc truy vấn trước.

Cấm viết "đã chứng minh" từ báo cáo này. Không nâng nhãn bằng chứng.
"""


def _interpret_target(r, min_cov):
    cov = float(r["qcov"])
    stops = int(r["n_stop_supported"])
    fs = int(r["n_frameshift"])
    if cov < min_cov:
        base = "có mảnh, độ phủ thấp (<{c:.0%})".format(c=min_cov)
    else:
        base = "có mảnh, độ phủ đạt"
    if stops or fs:
        return base + "; {s} codon dừng có chỗ dựa + {f} lệch khung ứng viên → kiểm thủ công, rồi đối chiếu Meredith 2014".format(s=stops, f=fs)
    return base + "; không thấy tổn thương trong vùng căn được → chưa kết luận gene nguyên vẹn (vùng ngoài HSP chưa xét)"


def _interpret_control(r, min_cov):
    cov = float(r["qcov"])
    stops = int(r["n_stop_supported"])
    fs = int(r["n_frameshift"])
    ok = cov >= min_cov and stops == 0 and fs == 0
    return ok, ("ĐẠT" if ok else "KHÔNG ĐẠT (cov={c:.2f}, stop={s}, fs={f})".format(c=cov, s=stops, f=fs))


def cmd_report(args):
    _, loci = _read_tsv(args.loci)
    _, lesions = _read_tsv(args.lesions)
    manifest_rows = _read_tsv(args.manifest)[1] if args.manifest else []
    today = args.date or _dt.date.today().isoformat()

    genes_expected = {}
    for m in manifest_rows:
        genes_expected[m["gene"]] = m
    for r in loci:
        genes_expected.setdefault(r["gene"], {"gene": r["gene"], "role": r["role"], "accession": r.get("accession", ""), "status": "ok"})

    best = {}
    for r in loci:
        g = r["gene"]
        if g not in best or float(r["sum_bitscore"]) > float(best[g]["sum_bitscore"]):
            best[g] = r
    n_loci_by_gene = {}
    for r in loci:
        n_loci_by_gene[r["gene"]] = n_loci_by_gene.get(r["gene"], 0) + 1

    out = []
    out.append("# B5 — quét mảnh vỡ gene răng (và GULO) trên bộ gene chào mào")
    out.append("")
    out.append("Ngày tạo: {d} · công cụ: `pipeline/p07_tooth_genes/tooth_scan.py report` · đầu vào: `{l}`, `{s}`".format(d=today, l=args.loci, s=args.lesions))
    out.append("")
    out.append("Chữ ký mọi kết quả: **STRONG_INFERENCE tối đa** — xem mục cuối trước khi trích dẫn bất kỳ dòng nào.")
    out.append("")

    # Doi chung
    out.append("## 1. Đối chứng dương")
    out.append("")
    controls = [g for g, m in genes_expected.items() if m.get("role") == ROLE_CONTROL]
    all_ok = True
    if not controls:
        all_ok = False
        out.append("**Không có gene đối chứng nào** trong lần chạy này → kết quả gene mục tiêu KHÔNG diễn giải được. Thêm `--controls` (vd MEPE, SPARC) rồi chạy lại.")
    else:
        out.append("| Gene | Accession | Locus tốt nhất | Độ phủ | Stop (chỗ dựa) | Lệch khung | Kết quả |")
        out.append("|---|---|---|---|---|---|---|")
        for g in sorted(controls):
            m = genes_expected[g]
            r = best.get(g)
            if r is None:
                all_ok = False
                out.append("| {g} | {a} | — | 0 | — | — | KHÔNG ĐẠT (không tìm thấy) |".format(g=g, a=m.get("accession", "")))
                continue
            ok, txt = _interpret_control(r, args.min_control_cov)
            all_ok = all_ok and ok
            out.append("| {g} | {a} | {s}:{x}-{y}({st}) | {c:.2f} | {k} | {f} | {t} |".format(
                g=g, a=r.get("accession", ""), s=r["sseqid"], x=r["s_start"], y=r["s_end"], st=r["strand"],
                c=float(r["qcov"]), k=r["n_stop_supported"], f=r["n_frameshift"], t=txt))
        out.append("")
        out.append("**Kết luận đối chứng:** " + ("ĐẠT — quy trình phát hiện được gene nguyên vẹn; dòng gene mục tiêu bên dưới đọc được." if all_ok
                                                else "KHÔNG ĐẠT — dòng gene mục tiêu bên dưới **không diễn giải được** cho tới khi sửa."))
    out.append("")

    # Muc tieu
    out.append("## 2. Gene mục tiêu")
    out.append("")
    out.append("| Gene | Accession truy vấn | Số locus | Locus tốt nhất | Độ phủ | %id max | e-value | Stop thô / có chỗ dựa | Lệch khung | Ranh exon | Đọc tạm |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|")
    targets = [g for g, m in genes_expected.items() if m.get("role") != ROLE_CONTROL]
    for g in sorted(targets):
        m = genes_expected[g]
        r = best.get(g)
        if m.get("status", "ok") != "ok":
            out.append("| {g} | — | — | — | — | — | — | — | — | — | truy vấn không có ở ngoài nhóm ({s}) → chưa thử được |".format(g=g, s=m.get("status")))
            continue
        if r is None:
            out.append("| {g} | {a} | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |".format(g=g, a=m.get("accession", "")))
            continue
        out.append("| {g} | {a} | {n} | {s}:{x}-{y}({st}) | {c:.2f} | {p} | {e} | {k0} / {k1} | {f} | {i} | {t} |".format(
            g=g, a=r.get("accession", ""), n=n_loci_by_gene.get(g, 0), s=r["sseqid"], x=r["s_start"], y=r["s_end"], st=r["strand"],
            c=float(r["qcov"]), p=r["max_pident"], e=r["best_evalue"], k0=r["n_stop_raw"], k1=r["n_stop_supported"],
            f=r["n_frameshift"], i=r["n_intron_like"], t=_interpret_target(r, args.min_target_cov)))
    out.append("")
    if not all_ok:
        out.append("> Đối chứng dương không đạt: bảng trên **chỉ để gỡ lỗi**, không dùng cho bất kỳ phát biểu nào.")
        out.append("")

    # Ton thuong chi tiet
    out.append("## 3. Tổn thương ứng viên (chỉ dòng `supported=1`; đầy đủ trong `{s}`)".format(s=args.lesions))
    out.append("")
    sup = [l for l in lesions if l.get("supported") == "1"]
    if not sup:
        out.append("_Không có tổn thương nào qua lọc chỗ dựa._")
    else:
        out.append("| Gene | Scaffold | Mạch | Loại | Vị trí aa truy vấn | Toạ độ nt | Khung | Ghi chú | Khớp Meredith 2014? |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        for l in sup:
            out.append("| {g} | {s} | {st} | {t} | {q} | {n} | {f1}→{f2} | {note} | _điền tay_ |".format(
                g=l["gene"], s=l["sseqid"], st=l["strand"], t=l["type"], q=l["qpos"], n=l["s_nt"],
                f1=l["frame_from"], f2=l["frame_to"], note=l["note"]))
    out.append("")
    out.append("Cột cuối để trống có chủ đích: so với bảng phụ lục Meredith 2014 là bước **người làm**, không tự động.")
    out.append("")

    # Buoc tiep
    out.append("## 4. Bước tiếp theo (theo thứ tự)")
    out.append("")
    out.append("1. `tooth_scan.py extract --fna <genome.fna.gz> --loci {l} --pad 2000 --translate loci_translate.txt` rồi đọc 3 khung dịch quanh mỗi tổn thương.".format(l=args.loci))
    out.append("2. Lấy bảng phụ lục Meredith 2014 (DOI 10.1126/science.1254390) → đối chiếu vị trí từng tổn thương; chỉ dòng **khớp** mới ghi \"cùng kiểu hỏng\".")
    out.append("3. Chạy lại cùng truy vấn trên một bộ gene sẻ chất lượng cao (sẻ vằn / sẻ ngô lớn — docs/05 B6) để tách lỗi lắp ráp khỏi tổn thương chung.")
    out.append("4. Ghi kết quả vào `research/synthesis/` với nhãn; **không** sửa `research/raw/`.")
    out.append("")
    out.append(MANDATORY_WARNING)
    with smart_open_write(args.out) as f:
        f.write("\n".join(out))
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser():
    p = argparse.ArgumentParser(
        prog="tooth_scan.py",
        description="B5: quét mảnh vỡ gene răng/GULO trên bộ gene chào mào bằng tblastn (xem README.md).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch-queries", help="NCBI esearch+efetch -> queries.faa + manifest (+registry). CẦN MẠNG.")
    f.add_argument("--genes", default=",".join(DEFAULT_TARGET_GENES), help="gene mục tiêu, phẩy (mặc định 6 gene răng + GULO)")
    f.add_argument("--controls", default=",".join(DEFAULT_CONTROL_GENES), help="gene đối chứng dương, phẩy (mặc định MEPE,SPARC)")
    f.add_argument("--organism", default="Alligator mississippiensis", help="loài ngoài nhóm (mặc định cá sấu Mỹ; D06 docs/02)")
    f.add_argument("--email", required=True, help="email cho NCBI E-utilities (bắt buộc)")
    f.add_argument("--api-key", default=None, help="NCBI API key (tuỳ chọn, tăng hạn mức)")
    f.add_argument("--all-proteins", action="store_true", help="không giới hạn refseq[Filter]")
    f.add_argument("--retmax", type=int, default=20)
    f.add_argument("--out", required=True, help="queries.faa")
    f.add_argument("--manifest", required=True, help="queries_manifest.tsv")
    f.add_argument("--registry", default=None, help="data/registry.csv để ghi 1 dòng (sha256); bỏ trống = không ghi")
    f.add_argument("--date", default=None, help="ghi đè ngày (YYYY-MM-DD), dùng cho test")
    f.set_defaults(func=cmd_fetch_queries)

    q = sub.add_parser("parse", help="bảng tblastn -> loci.tsv + lesions.tsv (+hits.tsv)")
    q.add_argument("--tblastn", required=True, help='đầu ra tblastn -outfmt "%s"' % OUTFMT)
    q.add_argument("--manifest", default=None, help="queries_manifest.tsv (ánh xạ qseqid -> gene/role)")
    q.add_argument("--max-evalue", type=float, default=1e-5)
    q.add_argument("--max-locus-gap", type=int, default=100000, help="nt tối đa giữa 2 HSP cùng locus (mặc định 100 kb)")
    q.add_argument("--flank", type=int, default=10, help="số aa khớp mỗi bên để codon dừng có chỗ dựa")
    q.add_argument("--flank-min-ident", type=float, default=0.5)
    q.add_argument("--max-q-gap", type=int, default=20, help="aa tối đa giữa 2 HSP kề nhau trên truy vấn để xét lệch khung")
    q.add_argument("--max-fs-gap", type=int, default=40, help="nt tối đa (có thể âm) giữa 2 HSP để gọi lệch khung ứng viên")
    q.add_argument("--min-intron", type=int, default=40, help="nt từ đó khoảng cách subject được coi là ranh exon (intron)")
    q.add_argument("--out-loci", required=True)
    q.add_argument("--out-lesions", required=True)
    q.add_argument("--out-hits", default=None)
    q.set_defaults(func=cmd_parse)

    e = sub.add_parser("extract", help="trích vùng locus (+pad) từ genome .fna(.gz), tuỳ chọn 3 khung dịch")
    e.add_argument("--fna", required=True)
    e.add_argument("--loci", required=True)
    e.add_argument("--pad", type=int, default=2000)
    e.add_argument("--role", default=None, help="chỉ trích role này (target/control)")
    e.add_argument("--out", required=True, help="loci.fa (theo hướng đọc gene)")
    e.add_argument("--translate", default=None, help="file văn bản 3 khung dịch")
    e.set_defaults(func=cmd_extract)

    t = sub.add_parser("translate", help="FASTA nucleotide -> 3 khung dịch, đánh dấu '*'")
    t.add_argument("--fasta", required=True)
    t.add_argument("--out", default="-")
    t.set_defaults(func=cmd_translate)

    r = sub.add_parser("report", help="loci.tsv + lesions.tsv -> report.md (có cảnh báo bắt buộc)")
    r.add_argument("--loci", required=True)
    r.add_argument("--lesions", required=True)
    r.add_argument("--manifest", default=None)
    r.add_argument("--min-control-cov", type=float, default=0.7)
    r.add_argument("--min-target-cov", type=float, default=0.3)
    r.add_argument("--date", default=None)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_report)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
