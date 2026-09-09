#!/usr/bin/env python3
"""dilution_sim.py -- brief X7: mo phong "pha nuoc duong" (introgression) --
sau bao nhieu doi lai nguoc (backcross) hoac giao phoi noi bo (isolated) thi
dau vet mot su kien lai voi bo gene "goc" tro nen khong the phat hien.

Cau hoi khoa hoc (xem research/briefs/X7-mo-phong-pha-loang-introgression.md):
chu du an mo ta mot mo hinh "bo gene goc" tron voi "bo gene ban dia" roi loang
dan qua cac doi lai nguoc, giong pha nuoc duong, giong khoai lang Nhat qua
F1->F2->F3->F4. Day la mo hinh chuan cua introgression trong di truyen quan
the. Hai dai luong suy giam theo HAI CACH KHAC NHAU:

  - TI LE TO TIEN (ancestry fraction) giam theo cap so nhan duoi lai nguoc
    lien tiep (~init_fraction x 0.5^doi), nhung KHONG doi (chi dao dong do
    troi dat) neu quan the lai giao phoi noi bo thay vi lai nguoc.
  - CHIEU DAI DOAN TO TIEN (tract length) giam do tai to hop cat nho dan --
    o CA HAI che do, vi tai to hop xay ra moi doi bat ke ti le tong the co
    giam hay khong. Doan dai = su kien lai gan day; doan ngan = su kien lai
    xa xua. Day la "dong ho" phan tu, va cung la gioi han: khi doan qua ngan
    se lan vao nhieu, khong phuong phap nao con phat hien duoc.

Mo hinh mo phong (GIA DINH RO RANG cua script -- xem README.md muc "Mo hinh
va gia dinh"): MOI replicate x MOI nhiem sac the la 1 ca the LUONG BOI that
su, giu DUNG HAI ban sao (homolog A, homolog B), moi ban sao 1 danh sach
doan co nguon goc "goc". Doi 0 khoi tao XAC DINH kieu con lai F1 (hoac trung
gian neu --init-fraction != 0.5): homolog A la doan goc dai
min(1,2xinit_fraction)xchrom_length_cm, homolog B la phan con lai (0 neu
init_fraction<=0.5). MOI doi, MOI ca the sinh 1 giao tu bang giam phan (meiosis):
sinh so diem tai to hop theo Poisson (ky vong chrom_length_cm/100 giao cheo/
nhiem sac the), voi moi khoang giua 2 diem cheo ke tiep chon ngau nhien (dong
xu cong bang) MOT TRONG HAI ban sao CUA CHINH NO (homolog A hoac B) -- day
moi dung nghia "chon ngau nhien mot trong hai ban sao" cua brief. Giao tu do
ghep voi 1 giao tu THU HAI de tao the he ke tiep: (che do backcross) giao tu
thu hai luon la cua 1 ca the ban dia THUAN (ca 2 homolog rong, nen giao tu
luon rong -- khong can giam phan that vi ket qua chac chan rong); (che do
isolated) giao tu thu hai la san pham giam phan CUA MOT CA THE KHAC (mate)
rut ngau nhien tu chinh quan the o THE HE HIEN TAI (mate cung tu giam phan
DOC LAP tu 2 homolog cua chinh no). Quan trong: NEU dung mo hinh don gian hon
(1 danh sach doan/ca the, khong tach 2 homolog) thi voi quan the KHOI TAO
DONG NHAT (moi ca the giong het nhau o doi 0 -- dung de --generations 0 tra
ve dung --init-fraction), che do isolated se BI DONG BANG vinh vien (giam
phan giua 2 ban sao GIONG HET NHAU luon tra ve chinh no) -- day la LY DO
BAT BUOC phai mo hinh 2 homolog rieng, da tu phat hien va sua trong qua
trinh viet script nay (xem README.md muc "Mo hinh va gia dinh").

3 lenh con, chi dung Python 3.11 stdlib (argparse/math/random/statistics),
KHONG numpy/scipy, KHONG chay docker, KHONG tai du lieu mang:

  simulate  Markov don gian tren moi nhiem sac the -> tracts.tsv (tung doan
            theo doi/replicate/nhiem sac the) + sim_stats.tsv (thong ke tong
            hop theo doi).
  detect    sim_stats.tsv + 2 nguong NGUOI DUNG DAT (--min-fraction,
            --min-tract-cm) -> detect.tsv (doi nao con "phat hien duoc").
  report    sim_stats.tsv + detect.tsv -> report.md: bang theo doi + dien
            giai + CANH BAO BAT BUOC ve gioi han mo hinh.

QUAN TRONG: day la mo hinh MINH HOA (bo qua chon loc, thay doi kich thuoc
quan the, bien thien toc do tai to hop, va gia dinh ban do di truyen dong
nhat) -- cho BAC DO LON, KHONG PHAI con so dung de cong bo. Xem README.md.
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
import sys
from contextlib import contextmanager

# ---------------------------------------------------------------------------
# Tien ich dung chung (doc/ghi streaming, dinh dang o -- giong quy uoc
# pipeline/p04_anomaly/composition_scan.py; moi loi trong pipeline/ tu chua,
# khong import cheo module, dung quy uoc hien co cua du an)
# ---------------------------------------------------------------------------


@contextmanager
def smart_open_read(path):
    """Mo doc: '-' = stdin. newline='' de tu xu ly CRLF/LF, dong 'r' +
    rstrip('\\r\\n') o noi goi."""
    if path == "-":
        yield sys.stdin
    else:
        f = open(path, "r", newline="", encoding="utf-8")
        try:
            yield f
        finally:
            f.close()


@contextmanager
def smart_open_write(path):
    """Mo ghi: '-' = stdout. newline='' de ghi LF thuan tuy tren moi he dieu
    hanh."""
    if path == "-":
        yield sys.stdout
    else:
        f = open(path, "w", newline="", encoding="utf-8")
        try:
            yield f
        finally:
            f.close()


def _fmt_cell(value):
    """Dinh dang 1 gia tri cho dong TSV: float -> 6 chu so co nghia, con lai
    -> str() thuong (giu nguyen so nguyen, khong them .0)."""
    if isinstance(value, float):
        return "{:.6g}".format(value)
    return str(value)


def _read_tsv(path):
    """Doc 1 file TSV -> (header_list, [(lineno, fields_list), ...]). Bo qua
    dong rong. lineno bat dau tu 2 (dong 1 la header) de bao loi ro rang."""
    with smart_open_read(path) as f:
        header_line = f.readline()
        if not header_line:
            return [], []
        header = header_line.rstrip("\r\n").split("\t")
        rows = []
        for lineno, raw_line in enumerate(f, 2):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            rows.append((lineno, line.split("\t")))
        return header, rows


def _require_columns(header, required, err_prefix):
    idx = {name: i for i, name in enumerate(header)}
    missing = [c for c in required if c not in idx]
    if missing:
        raise SystemExit(
            "{}: thieu cot {} (dung dung file do lenh phu hop trong "
            "dilution_sim.py xuat ra chua?)".format(err_prefix, missing)
        )
    return idx


# ---------------------------------------------------------------------------
# Header co dinh cho 3 file TSV
# ---------------------------------------------------------------------------

TRACTS_HEADER = ["replicate", "generation", "chrom", "start_cm", "end_cm", "len_cm"]
STATS_HEADER = [
    "generation",
    "mean_ancestry_fraction",
    "sd_ancestry_fraction",
    "n_tracts",
    "mean_tract_cm",
    "median_tract_cm",
    "max_tract_cm",
]
DETECT_HEADER = ["generation", "fraction_ok", "tract_ok", "detectable"]


def _read_stats_rows(path, err_prefix):
    header, raw_rows = _read_tsv(path)
    if not header:
        raise SystemExit("{}: rong (khong co dong header)".format(err_prefix))
    idx = _require_columns(header, STATS_HEADER, err_prefix)
    rows = []
    for lineno, fields in raw_rows:
        try:
            rows.append(
                {
                    "generation": int(fields[idx["generation"]]),
                    "mean_ancestry_fraction": float(fields[idx["mean_ancestry_fraction"]]),
                    "sd_ancestry_fraction": float(fields[idx["sd_ancestry_fraction"]]),
                    "n_tracts": int(fields[idx["n_tracts"]]),
                    "mean_tract_cm": float(fields[idx["mean_tract_cm"]]),
                    "median_tract_cm": float(fields[idx["median_tract_cm"]]),
                    "max_tract_cm": float(fields[idx["max_tract_cm"]]),
                }
            )
        except (ValueError, IndexError) as e:
            raise SystemExit(
                "{}: dong {} khong doc duoc: {!r} ({})".format(err_prefix, lineno, fields, e)
            )
    return rows


def _read_detect_rows(path, err_prefix):
    header, raw_rows = _read_tsv(path)
    if not header:
        raise SystemExit("{}: rong (khong co dong header)".format(err_prefix))
    idx = _require_columns(header, DETECT_HEADER, err_prefix)
    rows = []
    for lineno, fields in raw_rows:
        try:
            rows.append(
                {
                    "generation": int(fields[idx["generation"]]),
                    "fraction_ok": int(fields[idx["fraction_ok"]]),
                    "tract_ok": int(fields[idx["tract_ok"]]),
                    "detectable": int(fields[idx["detectable"]]),
                }
            )
        except (ValueError, IndexError) as e:
            raise SystemExit(
                "{}: dong {} khong doc duoc: {!r} ({})".format(err_prefix, lineno, fields, e)
            )
    return rows


# ---------------------------------------------------------------------------
# Loi Markov: Poisson (Knuth, khong numpy) + doan (tract) clip/merge + giao tu
# ---------------------------------------------------------------------------


def _sample_poisson(rng, lam):
    """So diem tai to hop ngau nhien ~ Poisson(lam), thuat toan tich luy
    tich cac so ngau nhien deu (0,1) den khi vuot e^-lam -- thuat toan pho
    bien trong nhieu giao trinh mo phong (thuong duoc gan cho Knuth); KHONG
    tu tra cuu duoc DOI cu the trong phien nay vi brief cam tai du lieu
    mang, nen KHONG gan nguon ngoai -- chi dung vi day la thuat toan chuan/
    pho bien va da duoc kiem chung gian tiep qua test (ti le trung binh cua
    mo phong khop cong thuc ky vong toan hoc). O(lam) phep nhan/vong lap nen
    CHI phu hop lam nho (vai don vi -- dung voi chrom_length_cm o thang cM
    thong thuong cua nhiem sac the that); KHONG danh cho lam lon (hang tram
    tro len) vi se cham."""
    if lam <= 0.0:
        return 0
    limit = math.exp(-lam)
    k = 0
    p = 1.0
    while True:
        k += 1
        p *= rng.random()
        if p <= limit:
            return k - 1


def _clip_tracts(tracts, lo, hi):
    """Tra ve phan cua danh sach doan (da sap xep tang dan, khong chong lan)
    nam trong khoang [lo, hi)."""
    out = []
    for s, e in tracts:
        cs = s if s > lo else lo
        ce = e if e < hi else hi
        if ce > cs:
            out.append((cs, ce))
    return out


def _merge_tracts(pieces):
    """Sap xep va gop cac doan lien ke/chong lan thanh dang chuan (khong
    chong lan, khong con 2 doan lien ke ma khong gop) -- tranh dem trung
    doan khi 2 khoang lien tiep cung chon 1 nguon va vo tinh noi lai."""
    if not pieces:
        return []
    pieces = sorted(pieces)
    merged = [pieces[0]]
    for s, e in pieces[1:]:
        ls, le = merged[-1]
        if s <= le:
            if e > le:
                merged[-1] = (ls, e)
        else:
            merged.append((s, e))
    return merged


def _meiosis(rng, homolog_a, homolog_b, lam, chrom_length_cm):
    """Sinh 1 'giao tu' bang giam phan tu HAI ban sao CUA CHINH MOT ca the
    (homolog_a, homolog_b): sinh diem tai to hop ~Poisson(lam) tren
    [0, chrom_length_cm], voi moi khoang giua 2 diem cheo ke tiep tung dong
    xu cong bang chon homolog_a hoac homolog_b, cat doan 'goc' cua ban sao do
    vao dung khoang, roi gop lai. Tra ve danh sach doan (start, end) da
    chuan hoa (sap xep, khong chong lan/lien ke). Goi HAM NAY DOC LAP moi
    lan can 1 giao tu (ke ca khi cung 1 ca the duoc rut lam mate cho nhieu
    ca the khac trong cung 1 doi) -- giong 1 ca the that sinh nhieu giao tu
    doc lap qua nhieu lan giam phan, KHONG tai su dung lai 1 giao tu da sinh
    truoc do cho muc dich khac."""
    k = _sample_poisson(rng, lam)
    breakpoints = sorted(rng.uniform(0.0, chrom_length_cm) for _ in range(k))
    bounds = [0.0] + breakpoints + [chrom_length_cm]
    pieces = []
    for i in range(len(bounds) - 1):
        lo, hi = bounds[i], bounds[i + 1]
        if hi <= lo:
            continue
        src = homolog_a if rng.random() < 0.5 else homolog_b
        pieces.extend(_clip_tracts(src, lo, hi))
    return _merge_tracts(pieces)


def _init_individual(init_fraction, chrom_length_cm):
    """Trang thai luong boi khoi tao (doi 0) cho 1 ca the/nhiem sac the,
    XAC DINH (khong ngau nhien) sao cho ti le to tien trung binh dung bang
    init_fraction: kieu F1 khi init_fraction=0.5 (homolog A goc TOAN BO,
    homolog B ban dia TOAN BO); tong quat hoa tuyen tinh cho init_fraction
    bat ky trong [0,1] bang cach lap day homolog A truoc roi moi sang homolog
    B. Tra ve (homolog_a, homolog_b)."""
    L = chrom_length_cm
    if init_fraction <= 0.5:
        len_a = 2.0 * init_fraction * L
        homolog_a = [(0.0, len_a)] if len_a > 0.0 else []
        homolog_b = []
    else:
        homolog_a = [(0.0, L)]
        len_b = (2.0 * init_fraction - 1.0) * L
        homolog_b = [(0.0, len_b)] if len_b > 0.0 else []
    return homolog_a, homolog_b


def _pick_mate_index(rng, r, replicates):
    """Chon chi so 1 replicate KHAC r, deu ngau nhien trong quan the hien
    tai (che do isolated). replicates == 1 la truong hop suy bien (quan the
    1 ca the) -- khong con ai khac de ghep, tu ghep voi chinh no (tu phoi)
    nhu 1 gia dinh don gian hoa, ghi ro trong README."""
    if replicates == 1:
        return r
    m = rng.randrange(replicates - 1)
    if m >= r:
        m += 1
    return m


# ---------------------------------------------------------------------------
# simulate: Markov tren moi nhiem sac the -> tracts.tsv + sim_stats.tsv
# ---------------------------------------------------------------------------


def add_simulate_parser(sub):
    p = sub.add_parser(
        "simulate",
        help="Mo phong Markov pha loang gene 'goc' qua cac doi (backcross|isolated) -> tracts.tsv + sim_stats.tsv",
        description=(
            "Mo phong Markov don gian tren moi nhiem sac the: MOI replicate x "
            "MOI nhiem sac the giu 1 danh sach doan co nguon goc 'goc' (mo hinh "
            "haploid-lineage don gian hoa -- xem docstring dau file va README.md "
            "muc 'Mo hinh va gia dinh'). Sinh doi 0 (trang thai khoi tao xac dinh: "
            "1 doan goc dai init_fraction*chrom_length_cm tren moi nhiem sac the), "
            "roi voi moi doi tiep theo: sinh diem tai to hop ~Poisson(chrom_length_cm/100), "
            "moi khoang giua 2 diem cheo chon ngau nhien (dong xu cong bang) ban sao "
            "cua chinh no (self) hoac cua ban phoi (mate) -- che do 'backcross': mate "
            "luon la ca the ban dia THUAN (0 doan goc); che do 'isolated': mate la 1 "
            "replicate KHAC rut ngau nhien tu chinh quan the o the he hien tai. Ghi "
            "tracts.tsv (moi dong 1 doan) va sim_stats.tsv (moi dong 1 doi, thong ke "
            "tong hop toan quan the)."
        ),
    )
    p.add_argument(
        "--generations", type=int, required=True, help="So doi mo phong (>=0; 0 = chi doi khoi tao)"
    )
    p.add_argument(
        "--mode",
        choices=["backcross", "isolated"],
        required=True,
        help="backcross = moi doi ghep voi ca the ban dia thuan; isolated = ghep voi 1 ca the khac trong chinh quan the",
    )
    p.add_argument(
        "--init-fraction",
        dest="init_fraction",
        type=float,
        default=0.5,
        help="Ti le 'goc' luc khoi tao (doi 0), vd 0.5 = con lai F1 (mac dinh 0.5)",
    )
    p.add_argument(
        "--chrom-length-cm",
        dest="chrom_length_cm",
        type=float,
        default=100.0,
        help="Chieu dai moi nhiem sac the tinh bang cM (mac dinh 100.0)",
    )
    p.add_argument(
        "--n-chrom", dest="n_chrom", type=int, default=30, help="So nhiem sac the mo phong (mac dinh 30)"
    )
    p.add_argument(
        "--replicates", type=int, default=100, help="So ca the (dong lai doc lap) mo phong (mac dinh 100)"
    )
    p.add_argument("--seed", type=int, default=42, help="Seed cho random.Random (mac dinh 42)")
    p.add_argument("--out", required=True, help="tracts.tsv dau ra, hoac - de ghi ra stdout")
    p.add_argument("--stats", required=True, help="sim_stats.tsv dau ra, hoac - de ghi ra stdout")
    p.set_defaults(func=cmd_simulate)


def _write_generation(fout, fstats, generation, state, chrom_length_cm):
    """Ghi tat ca doan cua the he 'generation' vao fout (tracts.tsv) va 1
    dong thong ke tong hop vao fstats (sim_stats.tsv). replicate/chrom trong
    tracts.tsv danh so tu 1 (de doc, khong phai chi so mang tu 0). Moi ca the
    la luong boi (2 homolog) -- doan cua CA HAI homolog deu duoc ghi ra
    (khong co cot phan biet homolog trong tracts.tsv, dung dac ta brief), va
    mau so tinh ti le to tien nhan 2 (2 ban sao/nhiem sac the/ca the)."""
    replicates = len(state)
    n_chrom = len(state[0]) if replicates else 0

    replicate_fracs = []
    all_lens = []
    for r in range(replicates):
        total_donor = 0.0
        for c in range(n_chrom):
            for homolog in state[r][c]:
                for s, e in homolog:
                    length = e - s
                    all_lens.append(length)
                    total_donor += length
                    fout.write(
                        "\t".join(
                            _fmt_cell(v) for v in (r + 1, generation, c + 1, s, e, length)
                        )
                        + "\n"
                    )
        chrom_total_len = n_chrom * chrom_length_cm * 2.0  # 2 homolog/nhiem sac the
        replicate_fracs.append(total_donor / chrom_total_len if chrom_total_len > 0 else 0.0)

    n_tracts = len(all_lens)
    mean_frac = statistics.fmean(replicate_fracs) if replicate_fracs else 0.0
    sd_frac = statistics.pstdev(replicate_fracs) if len(replicate_fracs) > 1 else 0.0
    mean_len = statistics.fmean(all_lens) if all_lens else 0.0
    median_len = statistics.median(all_lens) if all_lens else 0.0
    max_len = max(all_lens) if all_lens else 0.0

    fstats.write(
        "\t".join(
            _fmt_cell(v)
            for v in (generation, mean_frac, sd_frac, n_tracts, mean_len, median_len, max_len)
        )
        + "\n"
    )


def cmd_simulate(args):
    if args.generations < 0:
        raise SystemExit("simulate: --generations phai >= 0")
    if args.replicates < 1:
        raise SystemExit("simulate: --replicates phai >= 1")
    if args.n_chrom < 1:
        raise SystemExit("simulate: --n-chrom phai >= 1")
    if args.chrom_length_cm <= 0.0:
        raise SystemExit("simulate: --chrom-length-cm phai > 0")
    if not (0.0 <= args.init_fraction <= 1.0):
        raise SystemExit("simulate: --init-fraction phai trong [0,1]")

    rng = random.Random(args.seed)
    lam = args.chrom_length_cm / 100.0
    replicates = args.replicates
    n_chrom = args.n_chrom
    chrom_length_cm = args.chrom_length_cm

    # Doi 0: MOI replicate/nhiem sac the deu khoi tao XAC DINH (khong ngau
    # nhien) kieu F1 (hoac trung gian) tu _init_individual -- nen ti le to
    # tien trung binh o doi 0 KHOP DUNG BANG --init-fraction cho MOI replicate,
    # khong chi xap xi thong ke (xem README + test_generation_zero_matches_init_fraction).
    # state[r][c] = (homolog_a, homolog_b) -- ca the LUONG BOI that su.
    init_pair = _init_individual(args.init_fraction, chrom_length_cm)
    state = [
        [(list(init_pair[0]), list(init_pair[1])) for _c in range(n_chrom)] for _r in range(replicates)
    ]

    with smart_open_write(args.out) as fout, smart_open_write(args.stats) as fstats:
        fout.write("\t".join(TRACTS_HEADER) + "\n")
        fstats.write("\t".join(STATS_HEADER) + "\n")
        _write_generation(fout, fstats, 0, state, chrom_length_cm)

        for g in range(1, args.generations + 1):
            new_state = [[None] * n_chrom for _r in range(replicates)]
            for c in range(n_chrom):
                for r in range(replicates):
                    homolog_a, homolog_b = state[r][c]
                    self_gamete = _meiosis(rng, homolog_a, homolog_b, lam, chrom_length_cm)
                    if args.mode == "backcross":
                        # Ca the ban dia THUAN: ca 2 homolog rong -> giam
                        # phan chac chan cho giao tu rong, bo qua de khoi
                        # ton loi goi/tieu thu RNG vo ich.
                        mate_gamete = []
                    else:
                        mate_r = _pick_mate_index(rng, r, replicates)
                        mate_a, mate_b = state[mate_r][c]
                        mate_gamete = _meiosis(rng, mate_a, mate_b, lam, chrom_length_cm)
                    new_state[r][c] = (self_gamete, mate_gamete)
            state = new_state
            _write_generation(fout, fstats, g, state, chrom_length_cm)

    return 0


# ---------------------------------------------------------------------------
# detect: sim_stats.tsv + 2 nguong nguoi dung dat -> detect.tsv
# ---------------------------------------------------------------------------


def add_detect_parser(sub):
    p = sub.add_parser(
        "detect",
        help="Danh gia doi nao con 'phat hien duoc' theo 2 nguong NGUOI DUNG DAT -> detect.tsv",
        description=(
            "Doc sim_stats.tsv (tu lenh 'simulate'). Voi moi doi: fraction_ok=1 neu "
            "mean_ancestry_fraction >= --min-fraction; tract_ok=1 neu max_tract_cm "
            ">= --min-tract-cm (ca hai deu dung nguong KEM THEO, >=, khong phai > -- "
            "quy uoc cua script vi bien lien tuc nen xac suat trung dung nguong gan "
            "bang 0). detectable=1 CHI KHI CA HAI fraction_ok VA tract_ok deu =1. "
            "QUAN TRONG: --min-fraction va --min-tract-cm la THAM SO NGUOI DUNG DAT, "
            "KHONG PHAI hang so tu nhien -- gia tri thuc te phai lay tu tai lieu "
            "phuong phap that. research/raw/R11-introgression-ghost-lineage.md (brief "
            "R11, DA XONG) khong tim duoc 1 nguong pho quat da kiem chung truc tiep -- "
            "chi co moc tham khao UNVERIFIED (tong hop Perplexity, chua tu mo tung bai "
            "goc): ti le pha tron ~1-5 phan tram co the phat hien o do sau thoi gian vua "
            "phai, cong suat giam manh khi duoi 1-2 phan tram; rieng kich ban lai nguoc "
            "LIEN TUC (dung voi che do 'backcross' cua script nay) khong tim duoc con so "
            "mo phong rieng. Cho toi khi co so DA KIEM CHUNG, MOI ket luan tu lenh 'detect' CHI mang "
            "tinh minh hoa (xem README.md muc 'Nguong phat hien')."
        ),
    )
    p.add_argument("--stats", required=True, help="sim_stats.tsv dau vao (tu lenh 'simulate')")
    p.add_argument(
        "--min-tract-cm",
        dest="min_tract_cm",
        type=float,
        required=True,
        help="Nguong chieu dai doan toi thieu (cM) de con 'phat hien duoc' -- THAM SO NGUOI DUNG DAT, khong phai hang so tu nhien",
    )
    p.add_argument(
        "--min-fraction",
        dest="min_fraction",
        type=float,
        required=True,
        help="Nguong ti le to tien toi thieu de con 'phat hien duoc' -- THAM SO NGUOI DUNG DAT, khong phai hang so tu nhien",
    )
    p.add_argument("--out", required=True, help="detect.tsv dau ra, hoac - de ghi ra stdout")
    p.set_defaults(func=cmd_detect)


def cmd_detect(args):
    rows = _read_stats_rows(args.stats, "detect: --stats")
    with smart_open_write(args.out) as fout:
        fout.write("\t".join(DETECT_HEADER) + "\n")
        for row in rows:
            fraction_ok = 1 if row["mean_ancestry_fraction"] >= args.min_fraction else 0
            tract_ok = 1 if row["max_tract_cm"] >= args.min_tract_cm else 0
            detectable = 1 if (fraction_ok and tract_ok) else 0
            fout.write(
                "\t".join(_fmt_cell(v) for v in (row["generation"], fraction_ok, tract_ok, detectable))
                + "\n"
            )
    return 0


# ---------------------------------------------------------------------------
# report: sim_stats.tsv + detect.tsv -> report.md (bang + dien giai + canh bao)
# ---------------------------------------------------------------------------

MANDATORY_WARNING = """## CANH BAO BAT BUOC -- doc truoc khi dung so trong bao cao nay cho bat ky ket luan nao

Mo hinh mo phong nay (`dilution_sim.py`) **BO QUA**:

1. **Chon loc** (selection) -- khong co doan to tien nao duoc uu tien giu
   lai hay dao thai; moi vi tri tren nhiem sac the trung tinh nhu nhau.
2. **Thay doi kich thuoc quan the** (population size change) -- `--replicates`
   co dinh suot mo phong, khong co that co (bottleneck) hay tang truong quan
   the thuc nhu quan the chim ngoai tu nhien.
3. **Bien thien toc do tai to hop doc nhiem sac the** (recombination rate
   heterogeneity) -- dung DUNG MOT ti le trung binh `chrom_length_cm/100`
   giao cheo/nhiem sac the/doi, khong co "diem nong" hay vung it tai to hop
   gan tam dong (centromere) nhu nhiem sac the that.
4. **Ban do di truyen dong nhat** (uniform genetic map) -- 1 cM tuong ung
   xac suat tai to hop nhu nhau o moi vi tri, khong khop ban do di truyen
   thuc cua bat ky loai chim nao (kem ca chao mao).

Vi vay ket qua trong bao cao nay chi cho **BAC DO LON** (order of magnitude)
ve so doi can de mot su kien lai tro nen kho phat hien -- **KHONG PHAI con
so dung de cong bo hay ket luan dinh luong** cho truong hop that. Nguong
`--min-fraction` / `--min-tract-cm` cua lenh `detect` la **THAM SO NGUOI
DUNG DAT**, khong phai hang so sinh hoc/tu nhien -- gia tri thuc te (bao
nhieu cM la "con phat hien duoc" voi phuong phap nao) phai lay tu tai lieu
phuong phap that. `research/raw/R11-introgression-ghost-lineage.md` (brief
R11, DA XONG) khong tim ra 1 nguong pho quat da kiem chung truc tiep -- chi
co moc tham khao UNVERIFIED (tong hop Perplexity, chua tu mo tung bai goc);
cho toi khi co so DA KIEM CHUNG, cot `detectable` chi mang tinh **minh hoa**.
"""

BACKCROSS_VS_ISOLATED = """## Backcross vs isolated -- khac nhau ra sao

Hai che do chi khac nhau o **ban phoi (mate)** duoc ghep voi giao tu tai to
hop moi doi -- moi thu khac (tai to hop, chon self/mate tung khoang) giong
het nhau:

- **backcross** (lai nguoc): ban phoi LUON la 1 ca the ban dia THUAN (0 doan
  goc). Vi vay moi doi, xac suat 1 vi tri con giu goc chi con **mot nua** so
  voi doi truoc (mat theo cap so nhan, xap xi `init_fraction x 0.5^doi`) --
  giong pha loang lien tuc bang nuoc sach, hoac khoai lang Nhat lai nguoc
  F1->F2->F3->F4. Ti le to tien **giam don dieu** qua cac doi.
- **isolated** (giao phoi noi bo): ban phoi la **1 replicate KHAC rut ngau
  nhien** tu chinh quan the dang mo phong (ban than no cung mang mot phan
  gene "goc", ti le trung binh xap xi ti le CA quan the). Ve ky vong toan
  hoc, ti le to tien trung binh cua CA quan the **khong doi** qua cac doi --
  chi dao dong do troi dat ngau nhien (kich thuoc mau huu han). Day la diem
  khac biet then chot so voi backcross, va la co so cua phep thu phan biet 2
  che do trong `tests/test_dilution_sim.py`.

**CA HAI che do** deu bi tai to hop **cat nho doan to tien moi doi** (du ti
le tong the co giam hay khong), nen **chieu dai doan trung binh giam don
dieu qua cac doi o ca hai che do** -- day la "dong ho" phan tu cua su kien
lai: doan dai = lai gan day, doan ngan = lai xa xua.
"""


def add_report_parser(sub):
    p = sub.add_parser(
        "report",
        help="Xuat report.md: bang theo doi + dien giai du lieu + CANH BAO BAT BUOC ve gioi han mo hinh",
        description=(
            "Doc sim_stats.tsv (tu 'simulate') va detect.tsv (tu 'detect', noi theo "
            "cot 'generation'), xuat report.md gom: CANH BAO BAT BUOC (mo hinh bo qua "
            "chon loc/kich thuoc quan the/bien thien tai to hop/gia dinh ban do dong "
            "nhat -- chi cho bac do lon), phan dien giai TU DU LIEU THAT (khong hardcode: "
            "ti le dau/cuoi, doi phat hien duoc xa nhat), giai thich khai niem backcross "
            "vs isolated khac nhau ra sao, va bang day du theo doi."
        ),
    )
    p.add_argument("--stats", required=True, help="sim_stats.tsv dau vao (tu lenh 'simulate')")
    p.add_argument("--detect", required=True, help="detect.tsv dau vao (tu lenh 'detect')")
    p.add_argument("--out", required=True, help="report.md dau ra, hoac - de ghi ra stdout")
    p.set_defaults(func=cmd_report)


def cmd_report(args):
    stats_rows = _read_stats_rows(args.stats, "report: --stats")
    detect_rows = _read_detect_rows(args.detect, "report: --detect")
    detect_by_gen = {row["generation"]: row for row in detect_rows}

    # Tu sap xep lai theo generation tang dan -- khong am tham tin thu tu co
    # san trong file (giong quy uoc 'report' cua composition_scan.py).
    stats_rows = sorted(stats_rows, key=lambda r: r["generation"])

    lines = []
    lines.append("# Bao cao mo phong pha loang gene 'goc' qua cac doi lai (X7 dilution_sim)")
    lines.append("")
    lines.append(
        "Nguon: `{}` + `{}` -- {} dong du lieu.".format(args.stats, args.detect, len(stats_rows))
    )
    lines.append("")
    lines.append(MANDATORY_WARNING)

    lines.append("## Dien giai (tu chinh du lieu trong file nay)")
    lines.append("")
    if not stats_rows:
        lines.append("(khong co dong du lieu nao trong `--stats`.)")
    else:
        gens = [r["generation"] for r in stats_rows]
        fracs = [r["mean_ancestry_fraction"] for r in stats_rows]
        g0, gN = gens[0], gens[-1]
        f0, fN = fracs[0], fracs[-1]

        if len(stats_rows) < 2:
            lines.append(
                "Chi co 1 diem du lieu (doi {}, ti le to tien trung binh = {}) -- chua du "
                "de danh gia xu huong qua nhieu doi.".format(g0, _fmt_cell(f0))
            )
        else:
            ratio = (fN / f0) if f0 > 0 else float("nan")
            lines.append(
                "Du lieu chay tu doi {} den doi {}. Ti le to tien trung binh di tu {} "
                "(doi {}) xuong {} (doi {}) -- ti so cuoi/dau = {}.".format(
                    g0, gN, _fmt_cell(f0), g0, _fmt_cell(fN), gN, _fmt_cell(ratio)
                )
            )
            if f0 > 0 and fN <= f0 * 0.5:
                lines.append(
                    "Xu huong nay **giam manh, don dieu** -- phu hop kieu 'lai nguoc' "
                    "(backcross): moi doi ti le mat xap xi mot nua."
                )
            elif f0 > 0 and fN >= f0 * 0.7:
                lines.append(
                    "Xu huong nay **khong giam he thong** (con lai gan muc ban dau, chi "
                    "dao dong) -- phu hop kieu 'giao phoi noi bo' (isolated)."
                )
            else:
                lines.append(
                    "Xu huong nay giam nhung khong ro rang thuoc kieu nao trong khoang "
                    "doi da chay -- can chay them doi de phan biet backcross/isolated "
                    "bang mat."
                )

        detectable_gens = [
            g for g in gens if g in detect_by_gen and detect_by_gen[g]["detectable"] == 1
        ]
        if detectable_gens:
            last_detectable = max(detectable_gens)
            lines.append(
                "Voi nguong da chon o lenh `detect` (xem CANH BAO BAT BUOC ben tren), su "
                "kien lai **con phat hien duoc toi da den doi {}** trong du lieu nay (doi "
                "`detectable=1` xa nhat); tu doi sau do trong du lieu nay khong con dong "
                "nao `detectable=1`.".format(last_detectable)
            )
        else:
            lines.append(
                "Voi nguong da chon o lenh `detect`, **khong co doi nao** trong du lieu "
                "nay dat `detectable=1` (hoac su kien lai da duoi nguong tu doi {}, hoac "
                "nguong dat qua chat -- day la tham so nguoi dung dat, xem canh bao "
                "tren).".format(g0)
            )
    lines.append("")

    lines.append(BACKCROSS_VS_ISOLATED)

    lines.append("## Bang theo doi")
    lines.append("")
    display_cols = STATS_HEADER + ["fraction_ok", "tract_ok", "detectable"]
    lines.append("| " + " | ".join(display_cols) + " |")
    lines.append("|" + "---|" * len(display_cols))
    for row in stats_rows:
        d = detect_by_gen.get(row["generation"])
        frac_ok = d["fraction_ok"] if d else ""
        tract_ok = d["tract_ok"] if d else ""
        detectable = d["detectable"] if d else ""
        cells = [row[c] for c in STATS_HEADER] + [frac_ok, tract_ok, detectable]
        lines.append("| " + " | ".join(_fmt_cell(v) for v in cells) + " |")
    lines.append("")

    with smart_open_write(args.out) as fout:
        fout.write("\n".join(lines) + "\n")

    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser():
    parser = argparse.ArgumentParser(
        prog="dilution_sim.py",
        description=(
            "Brief X7: mo phong Markov 'pha nuoc duong' (introgression) -- sau bao "
            "nhieu doi lai nguoc (backcross) hoac giao phoi noi bo (isolated) thi dau "
            "vet mot su kien lai voi bo gene 'goc' tro nen khong the phat hien. 3 lenh "
            "con: simulate, detect, report. Chi dung Python 3.11 stdlib. MO HINH MINH "
            "HOA -- bo qua chon loc/kich thuoc quan the/bien thien tai to hop, chi cho "
            "bac do lon, KHONG PHAI con so dung de cong bo. Xem README.md."
        ),
        epilog=(
            "Vi du:\n"
            "  simulate --generations 20 --mode backcross --seed 42 --out tracts.tsv --stats sim_stats.tsv\n"
            "  detect   --stats sim_stats.tsv --min-tract-cm 5 --min-fraction 0.01 --out detect.tsv\n"
            "  report   --stats sim_stats.tsv --detect detect.tsv --out report.md\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)
    add_simulate_parser(sub)
    add_detect_parser(sub)
    add_report_parser(sub)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
