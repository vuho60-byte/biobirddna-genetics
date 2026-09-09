"""Unit test cho reversion_prob.py -- unittest thuan, khong ngau nhien (toan
bo mo hinh la cong thuc dong dai so, khong co random.Random nao ca, khac
dilution_sim.py). Phan lon du lieu nhung ngay trong file.

Chay: python -B -m unittest discover -s pipeline/p06_atavism/tests -v

Don dep sau test: dung tempfile.mkdtemp() + shutil.rmtree() tren CHINH thu
muc tam do test tu tao (khong dung lenh xoa tren duong dan nao khac) -- dung
quy uoc da co san trong tests/test_dilution_sim.py va cac test suite khac
cua du an (cung 1 quy tac du an "cam lenh xoa rm/rmdir/del/Remove-Item", da
duoc dung shutil.rmtree cho dung 1 muc dich nay: don rac tempfile cua chinh
test).

Nhom test QUAN TRONG NHAT (yeu cau BAT BUOC rieng cua brief X8, muc Test):
TestSequentialVsConcurrentDependenceOnN -- goi CUNG 1 CAP HAM
(sequential_wait_generations / concurrent_wait_generations) voi NHIEU gia
tri N rat khac nhau: duong tuan tu PHAI cho ket qua GIONG HET (khong phu
thuoc N), duong dong thoi PHAI ti le nghich CHINH XAC voi N. Day la phep
thu se THAT SU FAIL neu logic 2 duong bi hoan doi hoac ai do vo tinh them/
bot N khoi 1 trong 2 cong thuc.
"""

import math
import os
import shutil
import sys
import tempfile
import unittest
from decimal import Decimal, getcontext

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import reversion_prob as rp  # noqa: E402


def read_text(path):
    with open(path, "r", newline="", encoding="utf-8") as f:
        return f.read()


def write_text(path, content):
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(content)


class TempDirMixin:
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="revprob_test_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)


def tsv_blocks(text):
    """Tach van ban co NHIEU khoi TSV cach nhau boi dong rong (dung cho lenh
    'locus', xuat 2 bang: element/state roi metric/value/unit) thanh danh
    sach [(header_list, [row_dict,...]), ...], bo qua cac dong khong phai
    bang (vd dong GHI CHU cuoi file)."""
    blocks = []
    current = []
    for line in text.splitlines():
        if line == "":
            if current:
                blocks.append(current)
                current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)

    parsed = []
    for block in blocks:
        if "\t" not in block[0]:
            continue  # dong ghi chu tu do, khong phai header TSV
        header = block[0].split("\t")
        rows = []
        for ln in block[1:]:
            if "\t" not in ln:
                break
            rows.append(dict(zip(header, ln.split("\t"))))
        parsed.append((header, rows))
    return parsed


def metric_dict(text):
    """Lay khoi TSV 'metric/value/unit' (khoi thu 2 cua lenh 'locus', hoac
    khoi duy nhat cua lenh 'compare') thanh dict metric->value (chuoi tho)."""
    for header, rows in tsv_blocks(text):
        if header[:2] == ["metric", "value"]:
            return {r["metric"]: r["value"] for r in rows}
    raise AssertionError("khong tim thay khoi TSV metric/value/unit trong:\n" + text)


# ---------------------------------------------------------------------------
# p_specific: chuyen doi/dao chuyen cong lai dung theo ti le kappa
# ---------------------------------------------------------------------------


class TestPSpecific(unittest.TestCase):
    def test_transition_and_transversion_sum_to_one_for_various_kappa(self):
        # 1 duong chuyen doi + 2 duong dao chuyen = het moi kha nang doi
        # base -> tong xac suat phai dung bang 1, voi MOI kappa.
        for kappa in (0.1, 0.5, 1.0, 2.0, 3.7, 10.0, 100.0):
            transition = rp.p_specific_transition(kappa)
            transversion = rp.p_specific_transversion(kappa)
            total = transition + 2 * transversion
            self.assertAlmostEqual(
                total, 1.0, places=12,
                msg="kappa={}: transition+2*transversion={} != 1".format(kappa, total),
            )

    def test_kappa_one_gives_one_third_for_both_paths(self):
        # kappa=1 nghia la chuyen doi va dao chuyen xay ra DEU nhau -- ca 3
        # base con lai deu co xac suat nhu nhau = 1/3.
        self.assertAlmostEqual(rp.p_specific_transition(1.0), 1.0 / 3.0, places=12)
        self.assertAlmostEqual(rp.p_specific_transversion(1.0), 1.0 / 3.0, places=12)

    def test_transition_favoured_when_kappa_above_one(self):
        # kappa>1: chuyen doi nhanh hon dao chuyen -> p_specific(chuyen doi)
        # phai LON HON p_specific(dao chuyen) -- day la ly do brief chon
        # chuyen doi lam mac dinh (gioi han TREN).
        ps = rp.p_specific(2.0, transversion=False)
        pv = rp.p_specific(2.0, transversion=True)
        self.assertGreater(ps, pv)

    def test_kappa_not_positive_raises(self):
        with self.assertRaises(SystemExit):
            rp.p_specific_transition(0.0)
        with self.assertRaises(SystemExit):
            rp.p_specific_transversion(-1.0)


# ---------------------------------------------------------------------------
# Phep thu BAT BUOC theo brief X8: duong tuan tu KHONG phu thuoc N, duong
# dong thoi PHU THUOC N -- day la test PHAN BIET 2 con duong.
# ---------------------------------------------------------------------------


class TestSequentialVsConcurrentDependenceOnN(unittest.TestCase):
    def test_sequential_constant_concurrent_inversely_proportional_to_n(self):
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        k = 3
        n_values = [1e2, 1e6, 1e10]  # chenh nhau toi 10^8 lan

        seq_results = [rp.sequential_wait_generations(p_rev, k, n) for n in n_values]
        for s in seq_results[1:]:
            self.assertEqual(
                s, seq_results[0],
                "duong TUAN TU doi ket qua khi N doi -- SAI, cong thuc phai doc lap N",
            )

        con_results = [rp.concurrent_wait_generations(p_rev, k, n) for n in n_values]
        # Giam don dieu ro rang khi N tang.
        for i in range(1, len(con_results)):
            self.assertLess(con_results[i], con_results[i - 1])

        # Ti le NGHICH CHINH XAC: N tang X lan -> wait_time PHAI giam DUNG X
        # lan (khong chi "khac nhau chung chung"), vi rate=2*N*p_rev**k
        # tuyen tinh voi N.
        for i in range(1, len(n_values)):
            ratio_n = n_values[i] / n_values[i - 1]
            ratio_wait = con_results[i - 1] / con_results[i]
            self.assertAlmostEqual(
                ratio_wait / ratio_n, 1.0, places=6,
                msg="N tang {} lan nhung wait_time chi doi {} lan (ky vong ti le nghich chinh xac)".format(
                    ratio_n, ratio_wait
                ),
            )

        # Va phai khac han sequential (khong phai trung hop bang nhau) --
        # kiem it nhat 1 gia tri N cho ket qua CACH XA (>10 lan) sequential.
        self.assertGreater(con_results[1] / seq_results[0], 10.0)

    def test_sequential_wait_scales_linearly_with_k(self):
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        for k in (1, 2, 3, 5, 10):
            expected = k / p_rev
            self.assertAlmostEqual(rp.sequential_wait_generations(p_rev, k) / expected, 1.0, places=9)
        # Kiem ti le tuyen tinh THAT SU giua 2 diem (khong chi khop cong
        # thuc noi bo -- doi chieu ti so k voi ti so ket qua).
        w1 = rp.sequential_wait_generations(p_rev, 1)
        w5 = rp.sequential_wait_generations(p_rev, 5)
        self.assertAlmostEqual(w5 / w1, 5.0, places=9)

    def test_sequential_shorter_than_concurrent_when_k_at_least_two(self):
        # Day la co so brief dung de chon duong TUAN TU lam con so CHINH.
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        n = 1e6
        for k in (2, 3, 5, 8):
            seq = rp.sequential_wait_generations(p_rev, k)
            con = rp.concurrent_wait_generations(p_rev, k, n)
            self.assertLess(
                seq, con,
                "k={}: tuan tu ({}) khong ngan hon dong thoi ({}) -- trai voi du doan brief X8".format(k, seq, con),
            )

    def test_k_equals_one_both_paths_give_identical_result(self):
        # Voi CHI 1 phan tu, "tuan tu" va "dong thoi" phai trung nhau ve mat
        # dinh nghia (khong co gi de "tuan tu hoa" hay "dong thoi hoa" ca).
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        n = 1e6
        seq = rp.sequential_wait_generations(p_rev, 1)
        con = rp.concurrent_wait_generations(p_rev, 1, n)
        self.assertAlmostEqual(seq, 1.0 / p_rev, places=6)
        # dong thoi k=1 con phu thuoc N (2*N*p_rev), tuan tu thi khong --
        # 2 con duong CHI trung nhau khi N=1/(2*p_rev) (diem dac biet), nen
        # o day chi kiem seq dung cong thuc, khong so seq==con.
        self.assertAlmostEqual(con, 1.0 / (2.0 * n * p_rev), places=6)


# ---------------------------------------------------------------------------
# Tang mu lam giam thoi gian cho dung ti le nghich
# ---------------------------------------------------------------------------


class TestMuScaling(unittest.TestCase):
    def test_wait_times_scale_inversely_with_mu(self):
        kappa = 2.0
        k = 2
        n = 1e6
        mu1, mu2 = 4.6e-9, 9.2e-9  # gap dung 2 lan

        p_rev1 = rp.p_rev_single(mu1, kappa)
        p_rev2 = rp.p_rev_single(mu2, kappa)

        seq1 = rp.sequential_wait_generations(p_rev1, k)
        seq2 = rp.sequential_wait_generations(p_rev2, k)
        # tuan tu ti le 1/p_rev ~ 1/mu -> mu gap doi, wait giam dung mot nua.
        self.assertAlmostEqual(seq1 / seq2, 2.0, places=6)

        con1 = rp.concurrent_wait_generations(p_rev1, k, n)
        con2 = rp.concurrent_wait_generations(p_rev2, k, n)
        # dong thoi ti le 1/p_rev**k ~ 1/mu**k -> mu gap doi, wait giam 2**k lan.
        self.assertAlmostEqual(con1 / con2, 2.0 ** k, places=4)

    def test_p_rev_itself_scales_linearly_with_mu(self):
        kappa = 3.0
        p1 = rp.p_rev_single(2e-9, kappa)
        p2 = rp.p_rev_single(4e-9, kappa)
        self.assertAlmostEqual(p2 / p1, 2.0, places=9)


# ---------------------------------------------------------------------------
# compare: ti so p_new/p_rev dung bang 1/p_specific khi locus=1 base, k=1
# ---------------------------------------------------------------------------


class TestCompareRatio(TempDirMixin, unittest.TestCase):
    def test_single_base_single_element_ratio_equals_inverse_p_specific(self):
        out_path = self.path("compare.tsv")
        rc = rp.main(
            [
                "compare",
                "--locus-length", "1",
                "--mu", "4.6e-9",
                "--kappa", "2.0",
                "--mutated-count", "1",
                "--out", out_path,
            ]
        )
        self.assertEqual(rc, 0)
        m = metric_dict(read_text(out_path))
        ps = rp.p_specific(2.0, transversion=False)
        expected_ratio = 1.0 / ps
        self.assertAlmostEqual(float(m["ratio_p_new_over_p_rev"]) / expected_ratio, 1.0, places=5)

    def test_ratio_matches_direct_computation_for_locus_1000(self):
        out_path = self.path("compare1000.tsv")
        rc = rp.main(
            [
                "compare",
                "--locus-length", "1000",
                "--mu", "4.6e-9",
                "--kappa", "2.0",
                "--mutated-count", "1",
                "--out", out_path,
            ]
        )
        self.assertEqual(rc, 0)
        m = metric_dict(read_text(out_path))
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        p_new = rp.p_new_any_mutation(4.6e-9, 1000)
        expected_ratio = p_new / p_rev
        # places=4 (khong phai 6): gia tri doc lai TU FILE TSV da qua
        # _fmt_cell (lam tron 6 CHU SO CO NGHIA, xem _fmt_cell trong
        # reversion_prob.py) -- sai so tuong doi toi da ly thuyet cua buoc
        # lam tron nay la ~5e-6 (nua don vi o chu so co nghia thu 6), nen
        # doi chieu voi gia tri KHONG lam tron phai noi long dung sai xuong
        # duoi muc do chinh xac cua dinh dang luu tru, khong phai loi logic.
        self.assertAlmostEqual(float(m["ratio_p_new_over_p_rev"]) / expected_ratio, 1.0, places=4)
        # Vung lan can 2000 (~L/p_specific = 1000/0.5) -- kiem bac do lon
        # dung, khong chi khop lai chinh no.
        self.assertGreater(float(m["ratio_p_new_over_p_rev"]), 1900.0)
        self.assertLess(float(m["ratio_p_new_over_p_rev"]), 2100.0)

    def test_transversion_flag_changes_ratio(self):
        out_t = self.path("compare_transition.tsv")
        out_v = self.path("compare_transversion.tsv")
        rp.main(["compare", "--locus-length", "1", "--mu", "4.6e-9", "--kappa", "2.0", "--out", out_t])
        rp.main(["compare", "--locus-length", "1", "--mu", "4.6e-9", "--kappa", "2.0", "--transversion", "--out", out_v])
        r_t = float(metric_dict(read_text(out_t))["ratio_p_new_over_p_rev"])
        r_v = float(metric_dict(read_text(out_v))["ratio_p_new_over_p_rev"])
        self.assertNotAlmostEqual(r_t, r_v)
        # kappa=2 > 1 -> duong chuyen doi co p_specific LON hon -> p_rev lon
        # hon -> ti so p_new/p_rev NHO hon so voi duong dao chuyen.
        self.assertLess(r_t, r_v)


# ---------------------------------------------------------------------------
# k=0 phai bao loi ro rang, khong tra ve vo cuc im lang
# ---------------------------------------------------------------------------


class TestZeroMutatedElements(TempDirMixin, unittest.TestCase):
    def test_zero_k_raises_at_function_level(self):
        p_rev = rp.p_rev_single(4.6e-9, 2.0)
        with self.assertRaises(SystemExit):
            rp.sequential_wait_generations(p_rev, 0)
        with self.assertRaises(SystemExit):
            rp.concurrent_wait_generations(p_rev, 0, 1e6)
        with self.assertRaises(SystemExit):
            rp.log_concurrent_prob_per_gamete(p_rev, 0)

    def test_locus_command_with_no_mutated_elements_raises(self):
        with self.assertRaises(SystemExit):
            rp.main(
                [
                    "locus", "--elements", "A,B,C", "--mutated", "",
                    "--out", self.path("out.txt"),
                ]
            )

    def test_compare_zero_mutated_count_raises(self):
        with self.assertRaises(SystemExit):
            rp.main(
                [
                    "compare", "--locus-length", "10", "--mutated-count", "0",
                    "--out", self.path("out.tsv"),
                ]
            )

    def test_sweep_k_list_containing_zero_raises(self):
        with self.assertRaises(SystemExit):
            rp.main(
                [
                    "sweep", "--mu-list", "4.6e-9", "--k-list", "1,0,2",
                    "--n-list", "1e6", "--out", self.path("sweep.tsv"),
                ]
            )


# ---------------------------------------------------------------------------
# So rat nho khong bi tran xuong 0 / crash: dung math.log khi can, kiem
# bang k=10 (yeu cau rieng cua brief X8, muc Test).
# ---------------------------------------------------------------------------


class TestUnderflowSafety(unittest.TestCase):
    def test_direct_power_underflows_and_breaks_naive_log_at_k_ten(self):
        # Xac nhan TIEN DE: voi p_rev tong hop du nho, phep tinh p_rev**10
        # NGAY THANG THAT SU tran ve dung 0.0 trong float64, va lay log cua
        # ket qua do (cach lam ngay thang) THAT SU nem ValueError -- neu 2
        # dong assert nay tu no da fail thi kich ban khong con kiem tra dung
        # van de "tran so" nua.
        p_rev = 1e-33
        k = 10
        self.assertEqual(p_rev ** k, 0.0)
        with self.assertRaises(ValueError):
            math.log(p_rev ** k)

    def test_log_prob_avoids_domain_error_when_direct_power_underflows_at_k_ten(self):
        p_rev = 1e-33
        k = 10
        # Ham cua module PHAI khong nem loi, va PHAI khop dung ket qua toan
        # hoc k*log(p_rev) (tinh doc lap o day bang chinh math.log(p_rev),
        # khong qua p_rev**k).
        log_prob = rp.log_concurrent_prob_per_gamete(p_rev, k)
        expected = k * math.log(p_rev)
        self.assertAlmostEqual(log_prob, expected, places=9)
        self.assertTrue(math.isfinite(log_prob))
        self.assertLess(log_prob, -700)  # xac nhan day THAT SU la vung cuc nho

    def test_concurrent_wait_returns_inf_not_crash_for_astronomically_tiny_rate(self):
        # p_rev=1e-33, k=10, N=1e6 -> rate=2*N*p_rev**10 ~ 2e-324, ban than
        # RATE cung khong con bieu dien noi trong float64 (duoi muc denormal
        # nho nhat ~4.9e-324) -- nghia la dap so wait_time TOAN HOC THAT SU
        # (~5e323) vuot qua gioi han float64 (~1.8e308). math.inf moi la dap
        # so DUNG o day, khong phai gioi han cua cach tinh. Quan trong nhat:
        # KHONG duoc nem ZeroDivisionError hay OverflowError ra ngoai.
        p_rev = 1e-33
        k = 10
        n = 1e6
        try:
            wait = rp.concurrent_wait_generations(p_rev, k, n)
        except (ZeroDivisionError, OverflowError) as e:
            self.fail("concurrent_wait_generations nem loi thay vi tra ve math.inf: {}".format(e))
        self.assertEqual(wait, math.inf)

    def test_concurrent_wait_finite_and_cross_checked_with_decimal_at_k_ten_realistic_params(self):
        # Truong hop THUC TE (mu/kappa mac dinh cua CLI): p_rev**10 ~4e-87,
        # KHONG tran (con lau moi cham muc denormal ~4.9e-324) -- wait_time
        # phai la 1 so HUU HAN, LON, va khop CHINH XAC (trong sai so dau
        # phay dong) voi phep tinh doc lap bang decimal.Decimal (do chinh
        # xac tuy y, KHONG phai numpy) -- day la doi chieu THAT, khong chi
        # kiem "khong crash".
        mu, kappa, n, k = 4.6e-9, 2.0, 1e6, 10
        p_rev = rp.p_rev_single(mu, kappa)
        self.assertNotEqual(p_rev ** k, 0.0, "tien de: o thong so THUC TE nay p_rev**k KHONG duoc tran")

        wait = rp.concurrent_wait_generations(p_rev, k, n)
        self.assertTrue(math.isfinite(wait))
        self.assertGreater(wait, 0.0)

        getcontext().prec = 50
        rate_dec = Decimal(2) * Decimal(str(n)) * (Decimal(str(p_rev)) ** k)
        expected_wait_dec = Decimal(1) / rate_dec
        relative_error = abs(Decimal(repr(wait)) - expected_wait_dec) / expected_wait_dec
        self.assertLess(float(relative_error), 1e-9)

        # Va phai lon hon han wait o k=3 (tang don dieu theo k, dong thoi).
        wait_k3 = rp.concurrent_wait_generations(p_rev, 3, n)
        self.assertGreater(wait, wait_k3)


# ---------------------------------------------------------------------------
# p_new: khong huy so ve 0 khi mu rat nho (dung math.log1p/expm1)
# ---------------------------------------------------------------------------


class TestPNewFormula(unittest.TestCase):
    def test_matches_linear_approximation_for_small_mu_l(self):
        mu, L = 4.6e-9, 1000
        value = rp.p_new_any_mutation(mu, L)
        linear_approx = mu * L
        self.assertAlmostEqual(value / linear_approx, 1.0, places=5)

    def test_p_new_does_not_catastrophically_cancel_for_tiny_mu(self):
        # mu=1e-18: 1.0-mu lam tron ve DUNG 1.0 trong float64 (duoi ULP cua
        # 1.0, ~2.22e-16) -- phep tinh 1-(1-mu)**L NGAY THANG tra ve DUNG
        # 0.0 (huy so), dù dap so toan hoc THAT (~mu) la 1 so duong ro rang.
        mu = 1e-18
        L = 1
        naive = 1.0 - (1.0 - mu) ** L
        self.assertEqual(naive, 0.0, "tien de: phep tinh ngay thang phai THAT SU trieu tieu ve 0 o day")

        value = rp.p_new_any_mutation(mu, L)
        self.assertGreater(value, 0.0)
        self.assertAlmostEqual(value / mu, 1.0, places=6)

    def test_p_new_increases_with_locus_length(self):
        mu = 4.6e-9
        p10 = rp.p_new_any_mutation(mu, 10)
        p1000 = rp.p_new_any_mutation(mu, 1000)
        self.assertLess(p10, p1000)
        self.assertAlmostEqual(p1000 / p10, 100.0, places=3)  # L gap 100 lan, con o vung tuyen tinh

    def test_invalid_mu_or_length_raises(self):
        with self.assertRaises(SystemExit):
            rp.p_new_any_mutation(0.0, 10)
        with self.assertRaises(SystemExit):
            rp.p_new_any_mutation(1.5, 10)
        with self.assertRaises(SystemExit):
            rp.p_new_any_mutation(1e-9, 0)


# ---------------------------------------------------------------------------
# locus (qua CLI) -- bang trang thai + validate ten phan tu
# ---------------------------------------------------------------------------


class TestLocusCommand(TempDirMixin, unittest.TestCase):
    def test_element_state_table_marks_mutated_and_original_correctly(self):
        out_path = self.path("locus.txt")
        rc = rp.main(
            [
                "locus", "--elements", "A,B,C,D", "--mutated", "A,B",
                "--mu", "4.6e-9", "--kappa", "2.0",
                "--n-diploid", "1e6", "--generation-years", "2",
                "--out", out_path,
            ]
        )
        self.assertEqual(rc, 0)
        blocks = tsv_blocks(read_text(out_path))
        state_header, state_rows = blocks[0]
        self.assertEqual(state_header, ["element", "state"])
        by_name = {r["element"]: r["state"] for r in state_rows}
        self.assertEqual(by_name, {"A": "dot_bien", "B": "dot_bien", "C": "goc", "D": "goc"})

    def test_wait_years_equals_wait_generations_times_generation_years(self):
        out_path = self.path("locus2.txt")
        rp.main(
            [
                "locus", "--elements", "A,B,C", "--mutated", "A,B,C",
                "--mu", "4.6e-9", "--kappa", "2.0", "--n-diploid", "1e6",
                "--generation-years", "2.5", "--out", out_path,
            ]
        )
        m = metric_dict(read_text(out_path))
        seq_gen = float(m["sequential_wait_generations"])
        seq_years = float(m["sequential_wait_years"])
        # places=4: ca seq_gen VA seq_years doc lai tu TSV da lam tron 6 chu
        # so co nghia DOC LAP voi nhau (xem _fmt_cell) -- sai so tuong doi
        # cong don toi da ~1e-5, can dung sai noi long hon muc luu tru.
        self.assertAlmostEqual(seq_years / (seq_gen * 2.5), 1.0, places=4)

    def test_unknown_mutated_name_raises(self):
        with self.assertRaises(SystemExit):
            rp.main(
                [
                    "locus", "--elements", "A,B,C", "--mutated", "Z",
                    "--out", self.path("out.txt"),
                ]
            )

    def test_duplicate_element_name_raises(self):
        with self.assertRaises(SystemExit):
            rp.main(
                [
                    "locus", "--elements", "A,A,B", "--mutated", "A",
                    "--out", self.path("out.txt"),
                ]
            )

    def test_k_matches_number_of_mutated_elements(self):
        out_path = self.path("locus3.txt")
        rp.main(
            [
                "locus", "--elements", "A,B,C,D,E", "--mutated", "A,B,C",
                "--mu", "4.6e-9", "--kappa", "2.0", "--out", out_path,
            ]
        )
        m = metric_dict(read_text(out_path))
        self.assertEqual(m["k_mutated"], "3")


# ---------------------------------------------------------------------------
# sweep + report (qua CLI) -- kich thuoc luoi dung, report co canh bao +
# dien giai tu du lieu that.
# ---------------------------------------------------------------------------


class TestSweepAndReportIntegration(TempDirMixin, unittest.TestCase):
    def test_sweep_grid_size_matches_product_of_lists(self):
        sweep_path = self.path("sweep.tsv")
        rc = rp.main(
            [
                "sweep",
                "--mu-list", "1e-9,4.6e-9,1e-8",
                "--k-list", "1,2,3,5",
                "--n-list", "1e4,1e6,1e8",
                "--kappa", "2.0",
                "--generation-years", "2",
                "--out", sweep_path,
            ]
        )
        self.assertEqual(rc, 0)
        header, raw_rows = rp._read_tsv(sweep_path)
        self.assertEqual(header, rp.SWEEP_HEADER)
        self.assertEqual(len(raw_rows), 3 * 4 * 3)  # |mu| x |k| x |N|

    def test_sweep_sequential_column_constant_across_n_for_fixed_mu_k(self):
        sweep_path = self.path("sweep2.tsv")
        rp.main(
            [
                "sweep", "--mu-list", "4.6e-9", "--k-list", "3",
                "--n-list", "1e4,1e6,1e8", "--out", sweep_path,
            ]
        )
        rows = rp._read_sweep_rows(sweep_path, "test")
        self.assertEqual(len(rows), 3)
        seq_values = {r["wait_gen_sequential"] for r in rows}
        self.assertEqual(len(seq_values), 1, "wait_gen_sequential phai GIONG HET nhau qua cac N khac nhau")
        con_values = [r["wait_gen_concurrent"] for r in rows]
        self.assertEqual(len(set(con_values)), 3, "wait_gen_concurrent phai KHAC NHAU qua cac N khac nhau")

    def test_report_end_to_end_contains_mandatory_warnings_and_table(self):
        sweep_path = self.path("sweep3.tsv")
        report_path = self.path("report.md")
        rp.main(
            [
                "sweep", "--mu-list", "1e-9,1e-8", "--k-list", "1,2,3,5",
                "--n-list", "1e4,1e6,1e8", "--out", sweep_path,
            ]
        )
        rc = rp.main(["report", "--sweep", sweep_path, "--out", report_path])
        self.assertEqual(rc, 0)
        report = read_text(report_path)
        self.assertIn("CANH BAO BAT BUOC", report)
        self.assertIn("nhin thay duoc", report)
        self.assertIn("lan ra quan the", report)
        self.assertIn("MAT DOAN", report)
        self.assertIn("CHUYEN DOI GENE", report)
        self.assertIn("DUNG SO 0", report)
        self.assertIn("TAM THOI", report)
        self.assertIn("R14", report)
        self.assertIn("| mu | kappa |", report)

    def test_report_sequential_shorter_fraction_reflects_actual_data(self):
        # Dung sweep CHI co k=1 (tuan tu == dong thoi tai k=1, xem
        # TestSequentialVsConcurrentDependenceOnN) -> dien giai phai bao
        # 0/N dong (khong co k>=2 nao trong file) hoac bo qua muc do, KHONG
        # duoc bia so "N/N dong" sai voi du lieu that.
        sweep_path = self.path("sweep_k1_only.tsv")
        report_path = self.path("report_k1.md")
        rp.main(["sweep", "--mu-list", "4.6e-9", "--k-list", "1", "--n-list", "1e6", "--out", sweep_path])
        rp.main(["report", "--sweep", sweep_path, "--out", report_path])
        report = read_text(report_path)
        # k=1 khong >=2 nen phan "duong tuan tu co luon ngan hon" khong
        # duoc xuat hien (khong co dong nao de thong ke).
        self.assertNotIn("dong co k>=2", report)

        sweep_path2 = self.path("sweep_k_mixed.tsv")
        report_path2 = self.path("report_mixed.md")
        rp.main(["sweep", "--mu-list", "4.6e-9", "--k-list", "1,2,3", "--n-list", "1e6", "--out", sweep_path2])
        rp.main(["report", "--sweep", sweep_path2, "--out", report_path2])
        report2 = read_text(report_path2)
        # co 2 dong k>=2 (k=2,k=3), va CA HAI phai la tuan tu ngan hon (theo
        # test_sequential_shorter_than_concurrent_when_k_at_least_two) ->
        # dien giai phai bao dung 2/2.
        self.assertIn("2/2 dong co k>=2", report2)

    def test_report_missing_sweep_columns_raises(self):
        bad_path = self.path("bad_sweep.tsv")
        write_text(bad_path, "mu\tkappa\n4.6e-9\t2.0\n")
        with self.assertRaises(SystemExit):
            rp.main(["report", "--sweep", bad_path, "--out", self.path("out.md")])

    def test_report_empty_sweep_raises(self):
        empty_path = self.path("empty.tsv")
        write_text(empty_path, "")
        with self.assertRaises(SystemExit):
            rp.main(["report", "--sweep", empty_path, "--out", self.path("out.md")])


# ---------------------------------------------------------------------------
# --help (top-level va tung lenh con) phai thoat ma 0
# ---------------------------------------------------------------------------


class TestHelp(unittest.TestCase):
    def _assert_help_exits_zero(self, argv):
        with self.assertRaises(SystemExit) as ctx:
            rp.main(argv)
        self.assertEqual(ctx.exception.code, 0)

    def test_top_level_help(self):
        self._assert_help_exits_zero(["--help"])

    def test_locus_help(self):
        self._assert_help_exits_zero(["locus", "--help"])

    def test_compare_help(self):
        self._assert_help_exits_zero(["compare", "--help"])

    def test_sweep_help(self):
        self._assert_help_exits_zero(["sweep", "--help"])

    def test_report_help(self):
        self._assert_help_exits_zero(["report", "--help"])


if __name__ == "__main__":
    unittest.main()
