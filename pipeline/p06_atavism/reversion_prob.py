#!/usr/bin/env python3
"""reversion_prob.py -- brief X8: xac suat mot "gene co" (phan tu bi dot bien
lam hong chuc nang) HIEN LAI qua dot bien diem tu quay nguoc dung base goc,
va so sanh voi xac suat xuat hien MOT dot bien moi bat ky trong cung locus.

Cau hoi khoa hoc (xem research/briefs/X8-xac-suat-hien-lai-gen-co.md): chu du
an hoi (1) xac suat mot gene co HIEN LAI la bao nhieu, va (2) co BANG xac
suat xuat hien mot dot bien MOI khong. Brief dong khung o muc PHAN TU: 1
locus = danh sach phan tu dat ten A, B, C, ... Z, moi phan tu co trang thai
'goc' hoac 'dot_bien'. Ban dau k phan tu DAU TIEN o trang thai 'dot_bien',
phan con lai 'goc' -- cau tra loi phu thuoc hoan toan vao k (BAO NHIEU phan
tu phai quay nguoc), khong phu thuoc danh tinh/thu tu cu the cua tung phan
tu (mo hinh coi cac phan tu dot_bien la DONG NHAT/hoan doi duoc ve xac suat,
dung chung 1 cap mu/kappa/lua chon duong -- don gian hoa, xem README).

Cong thuc goc (tu brief, KHONG tu suy dan):
  p_rev = mu * p_specific            (xac suat 1 phan tu CU THE quay nguoc
                                       DUNG base goc, 1 the he, 1 ban sao NST)
  p_specific = kappa/(kappa+2)       neu duong CHUYEN DOI (transition)
  p_specific = 1/(kappa+2)           neu 1 duong DAO CHUYEN (transversion) cu the
  (kappa = ti le toc do chuyen doi/dao chuyen; moi base co 1 duong chuyen doi
  + 2 duong dao chuyen, nen p_specific(chuyen doi) + 2*p_specific(dao chuyen) = 1)

Mac dinh dung duong CHUYEN DOI (p_specific lon hon khi kappa>1) vi day la
TRUONG HOP THUAN LOI NHAT cho su hien lai -- con so dua ra la GIOI HAN TREN
cua kha nang hien lai, dung --transversion de doi sang duong con lai.

Hai con duong quay nguoc (ca hai deu cai, chon qua lenh con):
  DONG THOI: ca k phan tu cung quay nguoc trong 1 giao tu, 1 the he. Xac suat
    moi giao tu = p_rev**k. So ca ky vong moi the he trong CA quan the =
    2*N*p_rev**k (N = so ca the luong boi) -- PHU THUOC N.
  TUAN TU: tung phan tu quay nguoc rieng le, THEO THU TU, moi buoc TRUNG TINH
    (gene da hong san, 1 phan tu quay nguoc mot minh khong cho loi the chon
    loc nao) va phai CO DINH trong quan the truoc khi buoc sau xay ra. Thoi
    gian cho ky vong MOI BUOC ~ 1/p_rev the he, KHONG PHU THUOC N -- day la
    ket qua CHUAN cua ly thuyet trung tinh (Kimura): so dot bien trung tinh
    moi xuat hien moi the he trong quan the = 2N*p_rev, xac suat co dinh cua
    MOI dot bien trung tinh = 1/(2N), nhan lai N TRIET TIEU DUNG, cho toc do
    co dinh = p_rev/the he bat ke N. **GIA DINH -- CHUA TU MO NGUON XAC MINH
    TRONG PHIEN VIET SCRIPT NAY** (brief cam tai du lieu mang; day la ket
    qua giao trinh pho bien nhung KHONG duoc gan nhan CONFIRMED/STRONG_INFERENCE
    khi chua tu mo 1 nguon that trong chinh phien lam viec nay -- xem
    README.md muc 'Gia dinh' va tasks/lessons.md "Tri nho cua Main cung
    KHONG phai bang chung"). Tong k buoc DOC LAP -> tong ky vong = TONG cac
    ky vong (dung voi moi phan phoi, khong can gia thiet gi them) = k/p_rev.
  Duong TUAN TU cho thoi gian NGAN HON duong DONG THOI khi k>=2 (vi p_rev<1
  nen p_rev**k giam nhanh hon k lan) -- nen la CON SO CHINH tra loi cau hoi 1.

So sanh bat buoc (cau hoi 2 cua chu du an) -- xac suat xuat hien MOT dot bien
MOI bat ky o dau do trong locus dai L base, 1 the he, 1 ban sao NST:
  p_new = 1 - (1-mu)**L  ~=  L*mu
Ti so p_new/p_rev cho biet dot bien moi de xay ra hon bao nhieu lan so voi
quay nguoc DUNG 1 phan tu cu the.

On dinh so hoc (yeu cau rieng cua brief, xem README muc 'So rat nho'):
  - p_new tinh bang -math.expm1(L*math.log1p(-mu)) thay vi 1-(1-mu)**L truc
    tiep -- tranh HUY SO (catastrophic cancellation) khi mu rat nho: 1-mu co
    the lam tron ve DUNG 1.0 trong dau phay dong, khien (1-mu)**L cung tra
    ve 1.0 va p_new tra ve DUNG 0.0 dù dap so that la 1 so duong nho (~L*mu).
  - Xac suat duong DONG THOI (p_rev**k) tinh qua khong gian LOG (k*log(p_rev)
    roi moi exp lai) thay vi luy thua truc tiep -- p_rev**k co the tran ve
    DUNG 0.0 khi k lon/p_rev nho (float64 khong bieu dien duoc so duong nho
    hon ~4.9e-324), khien buoc tinh 1/rate tiep theo nem ZeroDivisionError
    dù dap so wait_time that (huu han hoac vo cuc) van tinh duoc dung qua
    log. Kiem bang k=10 (xem tests/test_reversion_prob.py).

4 lenh con, chi dung Python 3.11 stdlib (argparse/math/sys), KHONG numpy/
scipy, KHONG chay docker, KHONG tai du lieu mang, KHONG import module khac
trong pipeline/ (moi thu muc tu chua, dung quy uoc hien co cua du an -- xem
pipeline/p05_admixture/dilution_sim.py):

  locus    1 to hop phan tu (--elements/--mutated) -> bang trang thai tung
           phan tu + xac suat/thoi gian cho CA HAI duong, tinh ca the he VA
           nam.
  compare  So p_new (dot bien moi bat ky trong locus dai L) voi p_rev (quay
           nguoc dung k phan tu cu the) -> ti so -- tra loi THANG cau hoi 2.
  sweep    Luoi tham so (mu x k x N) -> sweep.tsv, vi hang so that (mu/kappa
           chuan cho chao mao) con dang cho brief R14
           (research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md, CHUA co
           luc viet script nay) -- ket qua phai BEN VUNG truoc 1 khoang gia
           tri, khong phu thuoc 1 con so dinh.
  report   sweep.tsv -> report.md: bang gon + 3 CANH BAO BAT BUOC (xuat hien
           # nhin thay # lan ra quan the; mo hinh bo qua chon loc/mat doan/
           chuyen doi gene; phan tu XOA HAN thi xac suat quay nguoc = 0,
           khong phai so nho).

QUAN TRONG: hang so mac dinh mu=4.6e-9, kappa=2.0 la VI DU MINH HOA lay
NGUYEN VAN tu brief X8, CHUA PHAI so da xac nhan rieng cho chao mao -- xem
README.md va DEFAULT_MU/DEFAULT_KAPPA duoi day.
"""

from __future__ import annotations

import argparse
import math
import sys
from contextlib import contextmanager

# ---------------------------------------------------------------------------
# Hang so TAM THOI -- lay nguyen van tu vi du CLI trong brief X8
# (research/briefs/X8-xac-suat-hien-lai-gen-co.md: "--mu 4.6e-9 --kappa 2.0
# ... --generation-years 2"), CHUA PHAI hang so da xac nhan rieng cho chao
# mao. Cho toi khi co research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md
# (dang cho subagent khac lam, CHUA xong luc viet script nay), MOI ket qua
# dung mac dinh nay CHI la MINH HOA thang bac do lon -- xem README.md.
# ---------------------------------------------------------------------------

DEFAULT_MU = 4.6e-9
DEFAULT_KAPPA = 2.0
DEFAULT_N_DIPLOID = 1e6
DEFAULT_GENERATION_YEARS = 2.0


# ---------------------------------------------------------------------------
# Tien ich dung chung (doc/ghi streaming, dinh dang o -- giong quy uoc
# pipeline/p05_admixture/dilution_sim.py; moi module trong pipeline/ tu chua,
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
    -> str() thuong (giu nguyen so nguyen, khong them .0). float('inf') in
    ra 'inf' -- gia tri hop le, nghia la thoi gian cho vuot qua kha nang bieu
    dien cua float64 (xem concurrent_wait_generations)."""
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
            "reversion_prob.py xuat ra chua?)".format(err_prefix, missing)
        )
    return idx


def _parse_name_list(raw, err_prefix):
    names = [x.strip() for x in raw.split(",") if x.strip() != ""]
    if not names:
        raise SystemExit("{}: danh sach rong".format(err_prefix))
    return names


def _parse_float_list(raw, err_prefix):
    try:
        values = [float(x.strip()) for x in raw.split(",") if x.strip() != ""]
    except ValueError as e:
        raise SystemExit("{}: khong doc duoc so thuc trong danh sach {!r} ({})".format(err_prefix, raw, e))
    if not values:
        raise SystemExit("{}: danh sach rong".format(err_prefix))
    return values


def _parse_int_list(raw, err_prefix):
    try:
        values = [int(x.strip()) for x in raw.split(",") if x.strip() != ""]
    except ValueError as e:
        raise SystemExit("{}: khong doc duoc so nguyen trong danh sach {!r} ({})".format(err_prefix, raw, e))
    if not values:
        raise SystemExit("{}: danh sach rong".format(err_prefix))
    return values


def _require_positive_int(k, name, err_prefix):
    """k=0 (khong phan tu nao dot bien) PHAI bao loi ro rang -- KHONG duoc
    tra ve vo cuc/0 im lang (yeu cau rieng cua brief X8, muc Test)."""
    if k < 1:
        raise SystemExit(
            "{}: {} phai >= 1 -- k=0 nghia la khong co phan tu nao can quay "
            "nguoc, khong dinh nghia duoc thoi gian cho".format(err_prefix, name)
        )


def _validate_elements_and_mutated(elements, mutated, err_prefix):
    if len(set(elements)) != len(elements):
        raise SystemExit("{}: --elements co ten trung lap: {}".format(err_prefix, elements))
    if len(set(mutated)) != len(mutated):
        raise SystemExit("{}: --mutated co ten trung lap: {}".format(err_prefix, mutated))
    elements_set = set(elements)
    unknown = [m for m in mutated if m not in elements_set]
    if unknown:
        raise SystemExit(
            "{}: --mutated co ten khong nam trong --elements: {}".format(err_prefix, unknown)
        )


# ---------------------------------------------------------------------------
# Loi mo hinh: p_specific, p_rev, hai con duong quay nguoc, p_new
# ---------------------------------------------------------------------------


def p_specific_transition(kappa):
    """p_specific tren duong CHUYEN DOI (transition): kappa/(kappa+2)."""
    if kappa <= 0.0:
        raise SystemExit("p_specific_transition: kappa phai > 0")
    return kappa / (kappa + 2.0)


def p_specific_transversion(kappa):
    """p_specific tren 1 duong DAO CHUYEN (transversion) cu the: 1/(kappa+2).
    Moi base co 2 duong dao chuyen -- ham nay tra ve xac suat cho DUNG 1
    trong 2 duong do, khong phai tong ca hai (xem docstring module: tong
    p_specific(chuyen doi) + 2*p_specific(dao chuyen) = 1, kiem trong test)."""
    if kappa <= 0.0:
        raise SystemExit("p_specific_transversion: kappa phai > 0")
    return 1.0 / (kappa + 2.0)


def p_specific(kappa, transversion=False):
    """Chon duong theo co --transversion. Mac dinh (transversion=False) dung
    duong CHUYEN DOI -- truong hop THUAN LOI NHAT (p_specific lon hon khi
    kappa>1), de con so dua ra la GIOI HAN TREN cua kha nang hien lai (yeu
    cau ro rang cua brief X8)."""
    return p_specific_transversion(kappa) if transversion else p_specific_transition(kappa)


def p_rev_single(mu, kappa, transversion=False):
    """p_rev = mu * p_specific -- xac suat 1 phan tu CU THE quay nguoc DUNG
    base goc, trong 1 the he, tren 1 ban sao nhiem sac the."""
    if not (0.0 < mu < 1.0):
        raise SystemExit("p_rev_single: mu phai nam trong khoang (0,1)")
    return mu * p_specific(kappa, transversion)


def sequential_wait_generations(p_rev, k, n_diploid=None):
    """Thoi gian cho ky vong (so the he) duong TUAN TU: k buoc DOC LAP, moi
    buoc thoi gian cho trung binh 1/p_rev the he (GIA DINH ket qua chuan ly
    thuyet trung tinh -- xem docstring dau file + README.md muc 'Gia dinh'),
    tong k buoc DOC LAP -> tong ky vong = TONG cac ky vong = k/p_rev (dung
    voi moi phan phoi thoi gian cho, KHONG can gia thiet phan phoi cu the --
    chi can moi buoc doc lap va co ky vong huu han).

    n_diploid NHAN VAO NHUNG KHONG DUNG trong cong thuc -- chi de giu CHUNG
    1 CHU KY GOI voi concurrent_wait_generations(p_rev, k, n_diploid), phuc
    vu (a) vong lap sweep goi dong loat ca 2 ham cho cung 1 luoi tham so, va
    (b) test PHAN BIET BAT BUOC theo brief X8
    (tests/test_reversion_prob.py::TestSequentialVsConcurrentDependenceOnN):
    goi CUNG 1 CAP HAM voi N khac nhau -- sequential PHAI cho ket qua GIONG
    HET (khong dung n_diploid), concurrent PHAI cho ket qua ti le nghich voi
    N. Neu ai do vo tinh them n_diploid vao cong thuc duoi day, test se THAT
    SU FAIL."""
    _require_positive_int(k, "k", "sequential_wait_generations")
    if p_rev <= 0.0:
        raise SystemExit("sequential_wait_generations: p_rev phai > 0")
    return k / p_rev


def log_concurrent_prob_per_gamete(p_rev, k):
    """log(p_rev**k) = k*log(p_rev) -- tinh TRUC TIEP tren log(p_rev) thay vi
    tinh p_rev**k roi moi lay log. Vi du p_rev=1e-33, k=10: p_rev**10 TRAN VE
    DUNG 0.0 trong float64 (nho hon so duong nho nhat bieu dien duoc,
    ~4.9e-324), khien math.log(p_rev**10) nem ValueError('math domain
    error'); ham nay khong bao gio gap loi do vi khong bao gio tinh p_rev**k
    truoc khi lay log (xem test_log_prob_avoids_domain_error_when_direct_power_underflows_at_k_ten)."""
    _require_positive_int(k, "k", "log_concurrent_prob_per_gamete")
    if p_rev <= 0.0:
        raise SystemExit("log_concurrent_prob_per_gamete: p_rev phai > 0")
    return k * math.log(p_rev)


def concurrent_prob_per_gamete(p_rev, k):
    """Xac suat CA k phan tu cung quay nguoc trong 1 giao tu, 1 the he =
    p_rev**k (gia dinh dong nhat: moi phan tu dot bien dung chung 1 cap
    mu/kappa/lua chon duong -- xem docstring dau file). Tinh qua
    exp(log_concurrent_prob_per_gamete(...)) -- vi 0<p_rev<1 nen log(p_rev)
    LUON AM, k*log(p_rev) LUON AM, exp() cua so am khong bao gio tran (chi co
    the ve 0.0 khi qua am, dung ve mat toan hoc: xac suat qua nho khong con
    phan biet duoc voi 0 o do chinh xac float64)."""
    return math.exp(log_concurrent_prob_per_gamete(p_rev, k))


def concurrent_expected_events_per_generation(p_rev, k, n_diploid):
    """So ca 'hien lai dong thoi' ky vong moi the he trong CA quan the =
    2*n_diploid*p_rev**k (2*N ban sao nhiem sac the trong quan the luong
    boi kich thuoc N -- dung cong thuc brief X8)."""
    if n_diploid <= 0.0:
        raise SystemExit("concurrent_expected_events_per_generation: n_diploid phai > 0")
    return 2.0 * n_diploid * concurrent_prob_per_gamete(p_rev, k)


def _exp_or_inf(x):
    """exp(x), tra ve math.inf thay vi nem OverflowError khi x qua lon
    (x > ~709.78, gioi han cua math.exp trong float64). Dung cho thoi gian
    cho khi rate qua nho: dap so toan hoc THAT SU vuot qua kha nang bieu dien
    cua float64 (~1.8e308) trong truong hop nay (xem
    concurrent_wait_generations) -- vo cuc LA dap so dung, khong phai gioi
    han cua cach tinh."""
    try:
        return math.exp(x)
    except OverflowError:
        return math.inf


def concurrent_wait_generations(p_rev, k, n_diploid):
    """Thoi gian cho ky vong (so the he) duong DONG THOI: xap xi 1/rate,
    rate = 2*n_diploid*p_rev**k = so ca ky vong xuat hien moi the he trong CA
    quan the (xem concurrent_expected_events_per_generation). Tinh trong
    khong gian LOG (log_rate = log(2*n_diploid) + k*log(p_rev), dung
    log_concurrent_prob_per_gamete) roi wait = exp(-log_rate) -- KHONG tinh
    p_rev**k truoc roi chia 1.0/rate truc tiep, vi voi k lon/p_rev nho,
    p_rev**k co the tran ve DUNG 0.0 (float64 khong bieu dien duoc so duong
    nho hon ~4.9e-324) khien phep chia 1.0/0.0 nem ZeroDivisionError, dù dap
    so wait_time THAT (huu han hoac vo cuc) van tinh duoc dung qua log.
    exp(-log_rate) tu no co the tran NGUOC LAI (qua LON, > ~1.8e308) khi rate
    qua nho gan muc p_rev**k khong con bieu dien noi nua -- luc do dap so
    toan hoc THAT SU vuot qua float64, tra ve math.inf la dung (xem
    _exp_or_inf va test_concurrent_wait_returns_inf_not_crash_for_astronomically_tiny_rate)."""
    _require_positive_int(k, "k", "concurrent_wait_generations")
    if p_rev <= 0.0:
        raise SystemExit("concurrent_wait_generations: p_rev phai > 0")
    if n_diploid <= 0.0:
        raise SystemExit("concurrent_wait_generations: n_diploid phai > 0")
    log_rate = math.log(2.0 * n_diploid) + log_concurrent_prob_per_gamete(p_rev, k)
    return _exp_or_inf(-log_rate)


def p_new_any_mutation(mu, locus_length):
    """Xac suat co IT NHAT 1 dot bien moi o dau do trong locus dai
    locus_length base, 1 the he, tren 1 ban sao nhiem sac the:
    p_new = 1 - (1-mu)**locus_length ~= locus_length*mu (khi locus_length*mu
    << 1). Tinh bang -math.expm1(locus_length*math.log1p(-mu)) thay vi
    1-(1-mu)**locus_length truc tiep -- tranh HUY SO (catastrophic
    cancellation) khi mu rat nho: 1-mu co the lam tron ve DUNG 1.0 trong dau
    phay dong (mu nho hon epsilon may ~2.22e-16), khien (1-mu)**locus_length
    cung tra ve DUNG 1.0 va p_new tra ve DUNG 0.0 dù dap so toan hoc that la
    1 so duong nho (~locus_length*mu) chu khong phai 0. math.log1p/math.expm1
    la ham CHUAN cua Python thu vien math, danh rieng cho truong hop nay (xem
    test_p_new_does_not_catastrophically_cancel_for_tiny_mu)."""
    if not (0.0 < mu < 1.0):
        raise SystemExit("p_new_any_mutation: mu phai nam trong khoang (0,1)")
    if locus_length < 1:
        raise SystemExit("p_new_any_mutation: locus_length phai >= 1")
    return -math.expm1(locus_length * math.log1p(-mu))


# ---------------------------------------------------------------------------
# locus: 1 to hop phan tu -> bang trang thai + xac suat/thoi gian cho ca 2 duong
# ---------------------------------------------------------------------------


def add_locus_parser(sub):
    p = sub.add_parser(
        "locus",
        help="1 locus (--elements/--mutated) -> bang trang thai + xac suat/thoi gian cho ca 2 duong quay nguoc",
        description=(
            "Nhan 1 danh sach ten phan tu (--elements) va tap con dang 'dot_bien' luc bat dau "
            "(--mutated, k = so phan tu). In bang trang thai tung phan tu, roi xac suat/thoi gian "
            "cho theo CA HAI duong quay nguoc (dong thoi -- phu thuoc N; tuan tu -- KHONG phu thuoc "
            "N, la CON SO CHINH khi k>=2, xem README), tinh ca the he VA nam."
        ),
    )
    p.add_argument("--elements", required=True, help="Danh sach ten phan tu, cach nhau dau phay, vd A,B,C")
    p.add_argument(
        "--mutated",
        required=True,
        help="Tap con cua --elements dang 'dot_bien' luc bat dau, cach nhau dau phay, vd A hoac A,B",
    )
    p.add_argument(
        "--mu", type=float, default=DEFAULT_MU,
        help="Toc do dot bien / vi tri / the he (mac dinh {} -- TAM THOI, cho R14)".format(DEFAULT_MU),
    )
    p.add_argument(
        "--kappa", type=float, default=DEFAULT_KAPPA,
        help="Ti le toc do chuyen doi/dao chuyen (mac dinh {} -- TAM THOI, cho R14)".format(DEFAULT_KAPPA),
    )
    p.add_argument(
        "--transversion", action="store_true",
        help="Dung duong DAO CHUYEN (p_specific thuong nho hon) thay vi duong CHUYEN DOI mac dinh (GIOI HAN TREN)",
    )
    p.add_argument(
        "--n-diploid", dest="n_diploid", type=float, default=DEFAULT_N_DIPLOID,
        help="Kich thuoc quan the luong boi N (mac dinh {:g})".format(DEFAULT_N_DIPLOID),
    )
    p.add_argument(
        "--generation-years", dest="generation_years", type=float, default=DEFAULT_GENERATION_YEARS,
        help="So nam moi the he (mac dinh {:g})".format(DEFAULT_GENERATION_YEARS),
    )
    p.add_argument("--out", default="-", help="File van ban dau ra, hoac - de ghi ra stdout (mac dinh -)")
    p.set_defaults(func=cmd_locus)


def cmd_locus(args):
    elements = _parse_name_list(args.elements, "locus: --elements")
    mutated = _parse_name_list(args.mutated, "locus: --mutated")
    _validate_elements_and_mutated(elements, mutated, "locus")
    k = len(mutated)
    _require_positive_int(k, "so phan tu --mutated", "locus")

    ps = p_specific(args.kappa, args.transversion)
    p_rev = p_rev_single(args.mu, args.kappa, args.transversion)

    seq_gen = sequential_wait_generations(p_rev, k, args.n_diploid)
    seq_years = seq_gen * args.generation_years
    con_gen = concurrent_wait_generations(p_rev, k, args.n_diploid)
    con_years = con_gen * args.generation_years
    con_prob = concurrent_prob_per_gamete(p_rev, k)
    con_expected = concurrent_expected_events_per_generation(p_rev, k, args.n_diploid)

    mutated_set = set(mutated)
    with smart_open_write(args.out) as f:
        f.write("\t".join(["element", "state"]) + "\n")
        for name in elements:
            state = "dot_bien" if name in mutated_set else "goc"
            f.write("{}\t{}\n".format(name, state))
        f.write("\n")

        f.write("\t".join(["metric", "value", "unit"]) + "\n")
        rows = [
            ("mu", args.mu, "per_site_per_generation"),
            ("kappa", args.kappa, "transition_transversion_ratio"),
            ("path_used", "transversion" if args.transversion else "transition", "category"),
            ("p_specific", ps, "probability"),
            ("p_rev", p_rev, "probability_per_generation_per_copy"),
            ("k_mutated", k, "count"),
            ("n_diploid", args.n_diploid, "individuals"),
            ("generation_years", args.generation_years, "years_per_generation"),
            ("sequential_wait_generations", seq_gen, "generations"),
            ("sequential_wait_years", seq_years, "years"),
            ("concurrent_prob_per_gamete", con_prob, "probability_per_gamete_per_generation"),
            ("concurrent_expected_events_per_generation", con_expected, "events_per_generation_population"),
            ("concurrent_wait_generations", con_gen, "generations"),
            ("concurrent_wait_years", con_years, "years"),
        ]
        for name, value, unit in rows:
            f.write("\t".join([name, _fmt_cell(value), unit]) + "\n")
        f.write("\n")
        f.write(
            "GHI CHU: sequential_wait_* la CON SO CHINH tra loi cau hoi 1 cua chu du an -- duong "
            "TUAN TU ngan hon duong DONG THOI khi k>=2 (xem README muc 'Hai con duong quay nguoc'). "
            "concurrent_* PHU THUOC n_diploid; sequential_* KHONG phu thuoc (GIA DINH, xem README "
            "muc 'Gia dinh').\n"
        )
    return 0


# ---------------------------------------------------------------------------
# compare: p_new (dot bien moi bat ky) vs p_rev (quay nguoc k phan tu cu the)
# ---------------------------------------------------------------------------


def add_compare_parser(sub):
    p = sub.add_parser(
        "compare",
        help="So p_new (dot bien moi bat ky trong locus dai L) voi p_rev (quay nguoc dung k phan tu) -> ti so",
        description=(
            "Tra loi THANG cau hoi 2 cua chu du an: xac suat xuat hien 1 dot bien MOI bat ky trong "
            "locus dai --locus-length base (p_new = 1-(1-mu)**L) co BANG xac suat quay nguoc DUNG "
            "--mutated-count phan tu cu the (p_rev = (mu*p_specific)**mutated_count, cong thuc DONG "
            "THOI voi k=mutated_count -- xem README) hay khong. In ti so p_new/p_rev."
        ),
    )
    p.add_argument("--locus-length", dest="locus_length", type=int, required=True, help="Do dai locus L, tinh bang base")
    p.add_argument(
        "--mu", type=float, default=DEFAULT_MU,
        help="Toc do dot bien / vi tri / the he (mac dinh {} -- TAM THOI, cho R14)".format(DEFAULT_MU),
    )
    p.add_argument(
        "--kappa", type=float, default=DEFAULT_KAPPA,
        help="Ti le toc do chuyen doi/dao chuyen (mac dinh {} -- TAM THOI, cho R14)".format(DEFAULT_KAPPA),
    )
    p.add_argument(
        "--transversion", action="store_true",
        help="Dung duong DAO CHUYEN thay vi duong CHUYEN DOI mac dinh (GIOI HAN TREN)",
    )
    p.add_argument(
        "--mutated-count", dest="mutated_count", type=int, default=1,
        help="k -- so phan tu dot bien dung de tinh p_rev (mac dinh 1)",
    )
    p.add_argument("--out", default="-", help="File van ban dau ra, hoac - de ghi ra stdout (mac dinh -)")
    p.set_defaults(func=cmd_compare)


def cmd_compare(args):
    _require_positive_int(args.mutated_count, "--mutated-count", "compare")
    ps = p_specific(args.kappa, args.transversion)
    p_rev_1 = p_rev_single(args.mu, args.kappa, args.transversion)
    p_rev_k = concurrent_prob_per_gamete(p_rev_1, args.mutated_count)
    p_new = p_new_any_mutation(args.mu, args.locus_length)

    if p_rev_k <= 0.0:
        raise SystemExit(
            "compare: p_rev^mutated_count qua nho, tran ve 0.0 trong float64 -- khong tinh duoc ti so "
            "p_new/p_rev (thu giam --mutated-count hoac tang --mu/--kappa)"
        )
    ratio = p_new / p_rev_k

    with smart_open_write(args.out) as f:
        f.write("\t".join(["metric", "value", "unit"]) + "\n")
        rows = [
            ("locus_length", args.locus_length, "bases"),
            ("mu", args.mu, "per_site_per_generation"),
            ("kappa", args.kappa, "transition_transversion_ratio"),
            ("path_used", "transversion" if args.transversion else "transition", "category"),
            ("p_specific", ps, "probability"),
            ("mutated_count", args.mutated_count, "count"),
            ("p_rev", p_rev_k, "probability_per_gamete_per_generation"),
            ("p_new", p_new, "probability_per_gamete_per_generation"),
            ("ratio_p_new_over_p_rev", ratio, "ratio"),
        ]
        for name, value, unit in rows:
            f.write("\t".join([name, _fmt_cell(value), unit]) + "\n")
    return 0


# ---------------------------------------------------------------------------
# sweep: luoi tham so (mu x k x N) -> sweep.tsv
# ---------------------------------------------------------------------------

SWEEP_HEADER = [
    "mu", "kappa", "transversion", "k", "n_diploid", "generation_years",
    "p_specific", "p_rev",
    "wait_gen_sequential", "wait_years_sequential",
    "wait_gen_concurrent", "wait_years_concurrent",
]


def add_sweep_parser(sub):
    p = sub.add_parser(
        "sweep",
        help="Luoi tham so (mu x k x N) -> sweep.tsv thoi gian cho ca 2 duong",
        description=(
            "Voi moi to hop (mu trong --mu-list, k trong --k-list, N trong --n-list): tinh "
            "p_specific/p_rev roi thoi gian cho ca 2 duong (the he + nam) -> 1 dong sweep.tsv. Can vi "
            "hang so mu/kappa THAT cho chao mao con dang cho brief R14 "
            "(research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md) -- ket qua phai BEN VUNG truoc 1 "
            "khoang gia tri, khong phu thuoc 1 con so dinh."
        ),
    )
    p.add_argument("--mu-list", dest="mu_list", required=True, help="Danh sach mu, cach nhau dau phay, vd 1e-9,4.6e-9,1e-8")
    p.add_argument(
        "--k-list", dest="k_list", default="1,2,3,5",
        help="Danh sach k (so phan tu dot bien), cach nhau dau phay (mac dinh 1,2,3,5)",
    )
    p.add_argument("--n-list", dest="n_list", required=True, help="Danh sach n_diploid, cach nhau dau phay, vd 1e4,1e6,1e8")
    p.add_argument(
        "--kappa", type=float, default=DEFAULT_KAPPA,
        help="Ti le toc do chuyen doi/dao chuyen (mac dinh {} -- TAM THOI, cho R14)".format(DEFAULT_KAPPA),
    )
    p.add_argument(
        "--transversion", action="store_true",
        help="Dung duong DAO CHUYEN thay vi duong CHUYEN DOI mac dinh (GIOI HAN TREN)",
    )
    p.add_argument(
        "--generation-years", dest="generation_years", type=float, default=DEFAULT_GENERATION_YEARS,
        help="So nam moi the he (mac dinh {:g})".format(DEFAULT_GENERATION_YEARS),
    )
    p.add_argument("--out", required=True, help="sweep.tsv dau ra, hoac - de ghi ra stdout")
    p.set_defaults(func=cmd_sweep)


def cmd_sweep(args):
    mu_list = _parse_float_list(args.mu_list, "sweep: --mu-list")
    k_list = _parse_int_list(args.k_list, "sweep: --k-list")
    n_list = _parse_float_list(args.n_list, "sweep: --n-list")
    for mu in mu_list:
        if not (0.0 < mu < 1.0):
            raise SystemExit("sweep: --mu-list co gia tri ngoai khoang (0,1): {}".format(mu))
    for k in k_list:
        _require_positive_int(k, "k (--k-list)", "sweep")
    for n in n_list:
        if n <= 0.0:
            raise SystemExit("sweep: --n-list co gia tri <= 0: {}".format(n))

    ps = p_specific(args.kappa, args.transversion)
    transversion_flag = 1 if args.transversion else 0

    with smart_open_write(args.out) as f:
        f.write("\t".join(SWEEP_HEADER) + "\n")
        for mu in mu_list:
            p_rev = p_rev_single(mu, args.kappa, args.transversion)
            for k in k_list:
                seq_gen = sequential_wait_generations(p_rev, k, None)
                seq_years = seq_gen * args.generation_years
                for n in n_list:
                    con_gen = concurrent_wait_generations(p_rev, k, n)
                    con_years = con_gen * args.generation_years
                    row = [
                        mu, args.kappa, transversion_flag, k, n, args.generation_years,
                        ps, p_rev, seq_gen, seq_years, con_gen, con_years,
                    ]
                    f.write("\t".join(_fmt_cell(v) for v in row) + "\n")
    return 0


# ---------------------------------------------------------------------------
# report: sweep.tsv -> report.md (bang gon + 3 CANH BAO BAT BUOC)
# ---------------------------------------------------------------------------

MANDATORY_WARNING = """## CANH BAO BAT BUOC -- doc truoc khi dung so trong bao cao nay cho bat ky ket luan nao

1. **"Xuat hien" khac "nhin thay duoc" khac "lan ra quan the".** Cong cu nay
   chi tinh xac suat MOT ban sao nhiem sac the quay lai DUNG base goc trong
   MOT the he. Neu tinh trang do la LAN (recessive), 1 ban quay nguoc o THE
   DI HOP (ket voi 1 ban con mang allele dot bien) se KHONG bieu hien ra
   kieu hinh -- can ca 2 ban cung mang allele goc (dong hop) moi thay duoc
   o 1 ca the, va can dot bien do LAN RA (tang tan so) trong quan the qua
   nhieu the he nua moi co y nghia hon 1 ca the don le. Ba cau hoi "xuat hien
   1 lan" / "nhin thay duoc o 1 ca the" / "lan ra quan the" la BA CAU HOI
   KHAC NHAU -- cong cu nay CHI tra loi cau dau tien.
2. **Mo hinh BO QUA**: chon loc (selection -- gia dinh moi phan tu trung
   tinh, khong loi the/bat loi khi quay nguoc), dot bien MAT DOAN (deletion
   -- neu phan tu bi xoa han thi KHONG con base nao de "quay nguoc" ve, xem
   canh bao 3 duoi day), va CHUYEN DOI GENE (gene conversion -- co che sao
   chep giua cac ban sao/paralog co the thay doi 1 base nhanh hon dot bien
   diem don thuan, khong duoc mo hinh hoa o day).
3. **Neu phan tu bi XOA HAN** (deletion that su, khong con base nao o vi tri
   do) thi xac suat quay nguoc qua CHINH CON DUONG dot bien diem (point
   mutation) nay bang **DUNG SO 0**, khong phai 1 so nho -- vi khong con gi
   de "dot bien tro lai dung base goc". Cong cu nay CHI mo hinh dot bien
   DIEM (thay 1 base bang 1 base khac o 1 vi tri VAN CON TON TAI); no khong
   ap dung cho truong hop mat doan/xoa gene.
"""


def _read_sweep_rows(path, err_prefix):
    header, raw_rows = _read_tsv(path)
    if not header:
        raise SystemExit("{}: rong (khong co dong header)".format(err_prefix))
    idx = _require_columns(header, SWEEP_HEADER, err_prefix)
    rows = []
    for lineno, fields in raw_rows:
        try:
            rows.append(
                {
                    "mu": float(fields[idx["mu"]]),
                    "kappa": float(fields[idx["kappa"]]),
                    "transversion": int(fields[idx["transversion"]]),
                    "k": int(fields[idx["k"]]),
                    "n_diploid": float(fields[idx["n_diploid"]]),
                    "generation_years": float(fields[idx["generation_years"]]),
                    "p_specific": float(fields[idx["p_specific"]]),
                    "p_rev": float(fields[idx["p_rev"]]),
                    "wait_gen_sequential": float(fields[idx["wait_gen_sequential"]]),
                    "wait_years_sequential": float(fields[idx["wait_years_sequential"]]),
                    "wait_gen_concurrent": float(fields[idx["wait_gen_concurrent"]]),
                    "wait_years_concurrent": float(fields[idx["wait_years_concurrent"]]),
                }
            )
        except (ValueError, IndexError) as e:
            raise SystemExit("{}: dong {} khong doc duoc: {!r} ({})".format(err_prefix, lineno, fields, e))
    return rows


def add_report_parser(sub):
    p = sub.add_parser(
        "report",
        help="sweep.tsv -> report.md: bang gon + 3 CANH BAO BAT BUOC ve gioi han mo hinh",
        description=(
            "Doc sweep.tsv (tu lenh 'sweep'), xuat report.md gom: 3 CANH BAO BAT BUOC (xuat hien vs "
            "nhin thay vs lan ra quan the; mo hinh bo qua chon loc/mat doan/chuyen doi gene; phan tu "
            "xoa han thi xac suat=0), 1 dong dien giai TU DU LIEU THAT (bao nhieu % dong k>=2 co tuan "
            "tu ngan hon dong thoi, khop du doan cua brief), va bang day du theo to hop (k, mu, N)."
        ),
    )
    p.add_argument("--sweep", required=True, help="sweep.tsv dau vao (tu lenh 'sweep')")
    p.add_argument("--out", required=True, help="report.md dau ra, hoac - de ghi ra stdout")
    p.set_defaults(func=cmd_report)


def cmd_report(args):
    rows = _read_sweep_rows(args.sweep, "report: --sweep")
    if not rows:
        raise SystemExit("report: --sweep khong co dong du lieu nao")

    rows_sorted = sorted(rows, key=lambda r: (r["k"], r["mu"], r["n_diploid"]))

    lines = []
    lines.append("# Bao cao xac suat gene co hien lai qua dot bien diem (X8 reversion_prob)")
    lines.append("")
    lines.append("Nguon: `{}` -- {} dong du lieu.".format(args.sweep, len(rows_sorted)))
    lines.append("")
    lines.append(MANDATORY_WARNING)

    lines.append("## Hang so dung trong bang duoi day la TAM THOI")
    lines.append("")
    lines.append(
        "`mu`/`kappa` trong bang la THAM SO NGUOI DUNG TRUYEN VAO khi chay lenh `sweep`, KHONG PHAI "
        "hang so sinh hoc da xac nhan cho chao mao. Gia tri mac dinh cua CLI (`mu={:g}`, `kappa={:g}`) "
        "lay NGUYEN VAN tu vi du minh hoa trong brief X8, cho TOI KHI "
        "`research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md` (dang cho, xem "
        "`research/briefs/R14-lai-to-va-hang-so-toc-do-dot-bien.md`) mang ve so da co nguon.".format(
            DEFAULT_MU, DEFAULT_KAPPA
        )
    )
    lines.append("")

    k2plus = [r for r in rows_sorted if r["k"] >= 2]
    if k2plus:
        n_seq_shorter = sum(1 for r in k2plus if r["wait_gen_sequential"] < r["wait_gen_concurrent"])
        lines.append("## Duong tuan tu co luon ngan hon duong dong thoi khi k>=2? (tu chinh du lieu nay)")
        lines.append("")
        lines.append(
            "{}/{} dong co k>=2 trong file nay co thoi gian cho TUAN TU ngan hon DONG THOI -- khop du "
            "doan cua brief X8 rang duong tuan tu la kich ban thuan loi nhat (con so CHINH) khi "
            "k>=2.".format(n_seq_shorter, len(k2plus))
        )
        lines.append("")

    lines.append("## Bang day du theo to hop (k, mu, n_diploid)")
    lines.append("")
    lines.append("| " + " | ".join(SWEEP_HEADER) + " |")
    lines.append("|" + "---|" * len(SWEEP_HEADER))
    for r in rows_sorted:
        cells = [r[c] for c in SWEEP_HEADER]
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
        prog="reversion_prob.py",
        description=(
            "Brief X8: xac suat mot 'gene co' (phan tu bi dot bien lam hong) HIEN LAI qua dot bien "
            "diem quay nguoc dung base goc, so voi xac suat xuat hien MOT dot bien moi bat ky trong "
            "locus. 4 lenh con: locus, compare, sweep, report. Chi Python 3.11 stdlib. Hang so mac "
            "dinh mu/kappa la TAM THOI, cho research/raw/R14-*.md. Xem README.md."
        ),
        epilog=(
            "Vi du:\n"
            "  locus   --elements A,B,C --mutated A --mu 4.6e-9 --kappa 2.0 --n-diploid 1e6 --generation-years 2\n"
            "  compare --locus-length 1000 --mu 4.6e-9 --kappa 2.0 --mutated-count 1\n"
            "  sweep   --mu-list 1e-9,4.6e-9,1e-8 --k-list 1,2,3,5 --n-list 1e4,1e6,1e8 --out sweep.tsv\n"
            "  report  --sweep sweep.tsv --out report.md\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)
    add_locus_parser(sub)
    add_compare_parser(sub)
    add_sweep_parser(sub)
    add_report_parser(sub)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
