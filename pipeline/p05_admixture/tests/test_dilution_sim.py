"""Unit test cho dilution_sim.py -- unittest thuan, seed co dinh, phan lon du
lieu nhung ngay trong file. Cac test mo hinh Markov that (simulate) dung
tham so nho/vua de nhanh nhung van du thong ke (>=200 replicate cho phep thu
cap so nhan, dung theo yeu cau brief X7).

Chay: python -B -m unittest discover -s pipeline/p05_admixture/tests -v

Don don sau test: dung tempfile.mkdtemp() + shutil.rmtree() tren CHINH thu
muc tam do test tu tao (khong dung lenh xoa tren duong dan nao khac) -- dung
quy uoc da co san trong tests/test_composition_scan.py, tests/test_short_elements.py,
tests/test_conserved_to_bulbul.py, tests/test_enrich.py (cung 1 quy tac du an
"cam lenh xoa rm/rmdir/del/Remove-Item", da duoc cac lan truoc dung
shutil.rmtree cho dung 1 muc dich nay: don rac tempfile cua chinh test).
"""

import math
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dilution_sim as ds  # noqa: E402


def read_text(path):
    with open(path, "r", newline="", encoding="utf-8") as f:
        return f.read()


def write_text(path, content):
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(content)


class TempDirMixin:
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="dilsim_test_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)


def tsv_rows(text):
    """Tach 1 khoi TSV (header + du lieu) thanh (header_list, [row_dict,...])."""
    lines = [ln for ln in text.splitlines() if ln != ""]
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        fields = ln.split("\t")
        rows.append(dict(zip(header, fields)))
    return header, rows


def standard_error(sd, n):
    return (sd / math.sqrt(n)) if n > 0 else 0.0


def normalize_ws(text):
    """Gop moi chuoi khoang trang/xuong dong lien tiep thanh 1 dau cach --
    dung de kiem tra 1 cum tu co xuat hien hay khong MA KHONG phu thuoc vao
    Markdown xuong dong mem (soft wrap) o giua cum tu do trong file .md (khi
    render Markdown that, xuong dong don trong 1 doan van duoc gap thanh dau
    cach, nen day la cach kiem dung voi noi dung DA RENDER, khong phai voi
    ky tu '\\n' tho trong file)."""
    return " ".join(text.split())


# ---------------------------------------------------------------------------
# _clip_tracts / _merge_tracts (don vi, khong qua CLI)
# ---------------------------------------------------------------------------


class TestClipAndMerge(unittest.TestCase):
    def test_clip_keeps_only_overlap(self):
        tracts = [(2.0, 5.0), (8.0, 10.0)]
        self.assertEqual(ds._clip_tracts(tracts, 0.0, 6.0), [(2.0, 5.0)])
        self.assertEqual(ds._clip_tracts(tracts, 3.0, 9.0), [(3.0, 5.0), (8.0, 9.0)])
        self.assertEqual(ds._clip_tracts(tracts, 10.0, 20.0), [])

    def test_clip_empty_source_returns_empty(self):
        self.assertEqual(ds._clip_tracts([], 0.0, 100.0), [])

    def test_merge_joins_adjacent_and_overlapping(self):
        pieces = [(5.0, 10.0), (10.0, 15.0), (20.0, 25.0)]
        self.assertEqual(ds._merge_tracts(pieces), [(5.0, 15.0), (20.0, 25.0)])

    def test_merge_overlapping_pieces(self):
        pieces = [(1.0, 5.0), (3.0, 8.0)]
        self.assertEqual(ds._merge_tracts(pieces), [(1.0, 8.0)])

    def test_merge_unsorted_input(self):
        pieces = [(20.0, 25.0), (1.0, 5.0)]
        self.assertEqual(ds._merge_tracts(pieces), [(1.0, 5.0), (20.0, 25.0)])

    def test_merge_empty(self):
        self.assertEqual(ds._merge_tracts([]), [])


# ---------------------------------------------------------------------------
# _sample_poisson (don vi, thong ke)
# ---------------------------------------------------------------------------


class TestSamplePoisson(unittest.TestCase):
    def test_zero_lambda_always_zero(self):
        rng = ds.random.Random(1)
        for _ in range(50):
            self.assertEqual(ds._sample_poisson(rng, 0.0), 0)

    def test_empirical_mean_close_to_lambda(self):
        # Kiem chung gian tiep thuat toan Knuth: trung binh mau ~5000 lan
        # rut phai gan lambda (dung sai so chuan cua trung binh Poisson =
        # sqrt(lambda/n)).
        rng = ds.random.Random(7)
        lam = 3.0
        n = 5000
        samples = [ds._sample_poisson(rng, lam) for _ in range(n)]
        mean = sum(samples) / n
        se = math.sqrt(lam / n)
        self.assertLess(abs(mean - lam), 6 * se, "trung binh mau {} qua xa lambda={}".format(mean, lam))

    def test_never_negative(self):
        rng = ds.random.Random(3)
        for _ in range(500):
            self.assertGreaterEqual(ds._sample_poisson(rng, 2.5), 0)


# ---------------------------------------------------------------------------
# _init_individual (don vi -- khoi tao doi 0, XAC DINH khong ngau nhien)
# ---------------------------------------------------------------------------


class TestInitIndividual(unittest.TestCase):
    def _avg_fraction(self, homolog_a, homolog_b, length):
        total = sum(e - s for s, e in homolog_a) + sum(e - s for s, e in homolog_b)
        return total / (2.0 * length)

    def test_half_init_fraction_is_classic_f1(self):
        # init_fraction=0.5 (mac dinh, con lai F1): 1 homolog goc TOAN BO, 1
        # homolog ban dia TOAN BO -- dung XAC DINH, khong ngau nhien.
        a, b = ds._init_individual(0.5, 100.0)
        self.assertEqual(a, [(0.0, 100.0)])
        self.assertEqual(b, [])
        self.assertAlmostEqual(self._avg_fraction(a, b, 100.0), 0.5, places=12)

    def test_low_init_fraction_only_fills_homolog_a(self):
        a, b = ds._init_individual(0.3, 100.0)
        self.assertEqual(b, [])
        self.assertAlmostEqual(self._avg_fraction(a, b, 100.0), 0.3, places=12)

    def test_high_init_fraction_fills_homolog_a_fully_then_b(self):
        a, b = ds._init_individual(0.7, 100.0)
        self.assertEqual(a, [(0.0, 100.0)])
        self.assertAlmostEqual(self._avg_fraction(a, b, 100.0), 0.7, places=12)

    def test_zero_and_one_edges(self):
        a0, b0 = ds._init_individual(0.0, 100.0)
        self.assertEqual((a0, b0), ([], []))
        a1, b1 = ds._init_individual(1.0, 100.0)
        self.assertAlmostEqual(self._avg_fraction(a1, b1, 100.0), 1.0, places=12)


# ---------------------------------------------------------------------------
# _meiosis (don vi) -- gom 1 test hoi quy TRUC TIEP cho loi da tim thay va
# sua trong qua trinh viet script: giam phan giua 2 ban sao GIONG HET NHAU
# phai tra ve DUNG chinh no (khong doi), nhung giam phan giua 1 ban sao day
# du + 1 ban sao rong phai cho ra KET QUA NGAU NHIEN (khac nhau giua cac lan
# rut).
# ---------------------------------------------------------------------------


class TestMeiosis(unittest.TestCase):
    def test_deterministic_given_seed(self):
        homolog_a = [(0.0, 100.0)]
        homolog_b = []
        r1 = ds._meiosis(ds.random.Random(99), homolog_a, homolog_b, 1.0, 100.0)
        r2 = ds._meiosis(ds.random.Random(99), homolog_a, homolog_b, 1.0, 100.0)
        self.assertEqual(r1, r2)

    def test_identical_homologs_reproduce_themselves_exactly(self):
        # HOI QUY: neu 2 'ban sao' giong het nhau, giam phan PHAI tra ve
        # dung chinh no bat ke diem cheo/dong xu ra sao (chon giua 2 thu
        # giong nhau khong the tao khac biet). Day chinh la nguyen nhan goc
        # cua loi "quan the isolated bi dong bang" da phat hien khi khoi tao
        # quan the dong nhat va dung mo hinh 1-ban-sao/ca-the: da sua bang
        # cach chuyen sang mo hinh luong boi that (2 homolog/ca the).
        same = [(10.0, 40.0), (60.0, 80.0)]
        for seed in (1, 2, 3, 4, 5):
            rng = ds.random.Random(seed)
            result = ds._meiosis(rng, same, list(same), 1.5, 100.0)
            self.assertEqual(result, same)

    def test_full_vs_empty_homolog_gives_roughly_half_length_on_average(self):
        # homolog_a day du (toan bo goc), homolog_b rong (toan bo ban dia):
        # ky vong chieu dai giao tu trung binh = chrom_length/2 (moi vi tri
        # 50% co hoi lay tu A). Kiem bang trung binh nhieu lan rut.
        homolog_a = [(0.0, 100.0)]
        homolog_b = []
        rng = ds.random.Random(11)
        n = 2000
        total_len = 0.0
        for _ in range(n):
            gamete = ds._meiosis(rng, homolog_a, homolog_b, 1.0, 100.0)
            total_len += sum(e - s for s, e in gamete)
        mean_len = total_len / n
        self.assertLess(abs(mean_len - 50.0), 3.0, "chieu dai giao tu trung binh {} qua xa 50".format(mean_len))

    def test_result_never_exceeds_chromosome_bounds(self):
        homolog_a = [(0.0, 100.0)]
        homolog_b = []
        rng = ds.random.Random(21)
        for _ in range(300):
            gamete = ds._meiosis(rng, homolog_a, homolog_b, 2.0, 100.0)
            for s, e in gamete:
                self.assertGreaterEqual(s, 0.0)
                self.assertLessEqual(e, 100.0)
                self.assertLess(s, e)


# ---------------------------------------------------------------------------
# _pick_mate_index (don vi)
# ---------------------------------------------------------------------------


class TestPickMateIndex(unittest.TestCase):
    def test_never_picks_self_when_population_over_one(self):
        rng = ds.random.Random(5)
        replicates = 10
        for r in range(replicates):
            for _ in range(50):
                m = ds._pick_mate_index(rng, r, replicates)
                self.assertNotEqual(m, r)
                self.assertGreaterEqual(m, 0)
                self.assertLess(m, replicates)

    def test_self_pairing_when_population_is_one(self):
        rng = ds.random.Random(5)
        self.assertEqual(ds._pick_mate_index(rng, 0, 1), 0)


# ---------------------------------------------------------------------------
# simulate (qua CLI, main()) -- day la nhom test QUAN TRONG NHAT: phai THAT
# SU phan biet duoc backcross vs isolated, khong chi chay khong loi.
# ---------------------------------------------------------------------------


class TestSimulateMathematicalProperties(TempDirMixin, unittest.TestCase):
    """Dung 1 cap tham so CHUNG (>=200 replicate theo dung yeu cau brief) cho
    ca backcross va isolated, roi kiem TAT CA tinh chat toan hoc brief doi
    hoi tu CUNG 2 lan chay nay (tranh mo phong lai nhieu lan, van du bao phu)."""

    GENERATIONS = 8
    REPLICATES = 200
    N_CHROM = 20
    CHROM_LENGTH_CM = 100.0
    SEED = 42
    INIT_FRACTION = 0.5

    def setUp(self):
        super().setUp()
        self.bc_tracts = self.path("bc_tracts.tsv")
        self.bc_stats = self.path("bc_stats.tsv")
        rc = ds.main(
            [
                "simulate",
                "--generations", str(self.GENERATIONS),
                "--mode", "backcross",
                "--init-fraction", str(self.INIT_FRACTION),
                "--chrom-length-cm", str(self.CHROM_LENGTH_CM),
                "--n-chrom", str(self.N_CHROM),
                "--replicates", str(self.REPLICATES),
                "--seed", str(self.SEED),
                "--out", self.bc_tracts,
                "--stats", self.bc_stats,
            ]
        )
        self.assertEqual(rc, 0)

        self.iso_tracts = self.path("iso_tracts.tsv")
        self.iso_stats = self.path("iso_stats.tsv")
        rc = ds.main(
            [
                "simulate",
                "--generations", str(self.GENERATIONS),
                "--mode", "isolated",
                "--init-fraction", str(self.INIT_FRACTION),
                "--chrom-length-cm", str(self.CHROM_LENGTH_CM),
                "--n-chrom", str(self.N_CHROM),
                "--replicates", str(self.REPLICATES),
                "--seed", str(self.SEED),
                "--out", self.iso_tracts,
                "--stats", self.iso_stats,
            ]
        )
        self.assertEqual(rc, 0)

        _, self.bc_rows = tsv_rows(read_text(self.bc_stats))
        _, self.iso_rows = tsv_rows(read_text(self.iso_stats))
        self.assertEqual(len(self.bc_rows), self.GENERATIONS + 1)
        self.assertEqual(len(self.iso_rows), self.GENERATIONS + 1)

    # -- Phep thu 1: backcross giam DON DIEU va khop cap so nhan trong sai
    # so thong ke (SE tinh tu chinh sd_ancestry_fraction/sqrt(replicates)
    # bao cao trong file, KHONG hardcode sai so tuyet doi). -------------------

    def test_backcross_fraction_decays_geometrically_and_monotonically(self):
        fracs = [float(r["mean_ancestry_fraction"]) for r in self.bc_rows]
        sds = [float(r["sd_ancestry_fraction"]) for r in self.bc_rows]

        # Don dieu: MOI doi phai <= doi truoc (khong duoc tang).
        for g in range(1, len(fracs)):
            self.assertLessEqual(
                fracs[g], fracs[g - 1],
                "backcross: ti le to tien TANG tu doi {} sang doi {} ({} -> {})".format(
                    g - 1, g, fracs[g - 1], fracs[g]
                ),
            )

        # Khop cong thuc ky vong init_fraction * 0.5**g trong sai so thong ke
        # (6 lan sai so chuan cua trung binh -- rat rong, van dam bao phat
        # hien duoc loi logic vi hieu ung mong doi la nhan doi/doi moi doi,
        # lon hon nhieu lan SE).
        for g, (frac, sd) in enumerate(zip(fracs, sds)):
            expected = self.INIT_FRACTION * (0.5 ** g)
            se = standard_error(sd, self.REPLICATES)
            tol = max(1e-9, 6 * se)
            self.assertLess(
                abs(frac - expected), tol,
                "backcross doi {}: quan sat {} qua xa ky vong {} (SE={}, tol={})".format(
                    g, frac, expected, se, tol
                ),
            )

    # -- Phep thu 2 (BAT BUOC theo brief -- phan biet 2 che do): isolated
    # KHONG giam he thong, chi dao dong; VA gia tri cuoi cua isolated PHAI
    # cach xa gia tri backcross du doan cung doi (cap so nhan), khong chi
    # "gan bang nhau tinh co". ------------------------------------------------

    def test_isolated_fraction_does_not_decay_and_differs_sharply_from_backcross(self):
        iso_fracs = [float(r["mean_ancestry_fraction"]) for r in self.iso_rows]

        # Khong giam he thong: MOI doi phai nam trong 1 bang RONG quanh
        # init_fraction (bien do thuc do duoc trong du lieu nay chi ~0.01 --
        # dung bang +/-0.2 la RAT RONG, chi de tranh flaky, van thua suc bat
        # duoc loi "isolated roi vao suy giam cap so nhan giong backcross").
        band = 0.2
        for g, frac in enumerate(iso_fracs):
            self.assertLess(
                abs(frac - self.INIT_FRACTION), band,
                "isolated doi {}: ti le {} lech qua xa init_fraction={} (>{})".format(
                    g, frac, self.INIT_FRACTION, band
                ),
            )

        # Phep thu PHAN BIET chinh: o doi cuoi, isolated phai CAO HON HAN
        # (>= 20 lan) so voi gia tri backcross se co o CUNG doi do -- day la
        # phep thu se THAT SU FAIL neu ai do vo tinh trao doi logic backcross
        # <-> isolated, hoac isolated bi lam giam ho theo cap so nhan.
        g_last = self.GENERATIONS
        backcross_predicted_at_last = self.INIT_FRACTION * (0.5 ** g_last)
        iso_last = iso_fracs[-1]
        self.assertGreater(
            iso_last,
            20 * backcross_predicted_at_last,
            "isolated doi {} = {} KHONG cao hon ro ret du doan backcross {} -- khong phan biet duoc 2 che do".format(
                g_last, iso_last, backcross_predicted_at_last
            ),
        )

    # -- Phep thu 3: chieu dai doan trung binh giam don dieu o CA HAI che do. --

    def test_mean_tract_length_decreases_monotonically_in_both_modes(self):
        for label, rows in (("backcross", self.bc_rows), ("isolated", self.iso_rows)):
            lens = [float(r["mean_tract_cm"]) for r in rows if int(r["n_tracts"]) > 0]
            self.assertGreaterEqual(len(lens), 2, "{}: qua it doi co doan de kiem xu huong".format(label))
            for g in range(1, len(lens)):
                self.assertLessEqual(
                    lens[g], lens[g - 1],
                    "{}: chieu dai doan trung binh TANG tu vi tri {} sang {} ({} -> {})".format(
                        label, g - 1, g, lens[g - 1], lens[g]
                    ),
                )

    # -- Phep thu 4: khong doan nao vuot [0, chrom_length_cm]; tong chieu dai
    # doan = mean_ancestry_fraction * replicates * n_chrom * chrom_length_cm
    # * 2 (luong boi: 2 homolog/nhiem sac the) trong sai so lam tron 6 chu so
    # co nghia cua dinh dang TSV. ---------------------------------------------

    def test_tracts_within_bounds_and_length_matches_fraction(self):
        for label, tracts_path, stats_rows in (
            ("backcross", self.bc_tracts, self.bc_rows),
            ("isolated", self.iso_tracts, self.iso_rows),
        ):
            _, trows = tsv_rows(read_text(tracts_path))
            sum_len_by_gen = {}
            for row in trows:
                s = float(row["start_cm"])
                e = float(row["end_cm"])
                length = float(row["len_cm"])
                self.assertGreaterEqual(s, 0.0, "{}: start_cm < 0".format(label))
                self.assertLessEqual(e, self.CHROM_LENGTH_CM, "{}: end_cm > chrom_length_cm".format(label))
                self.assertLess(s, e, "{}: doan rong/am".format(label))
                # start_cm/end_cm/len_cm deu duoc lam tron DOC LAP ve 6 chu so
                # CO NGHIA boi _fmt_cell truoc khi ghi TSV (xem _fmt_cell) --
                # sai so lam tron cua s/e ti le voi CHINH gia tri s/e (co the
                # lon hon nhieu so voi do dai doan e-s khi doan ngan nhung
                # nam xa vi tri 0, vd s=45.3 e=46.0: sai so lam tron cua s,e
                # ti le ~45 chu khong phai ~0.7) -- nen dung dung sai tuyet
                # doi tinh tu |s|+|e|, KHONG PHAI tu |e-s|, va KHONG dung
                # assertAlmostEqual(places=6) (don vi thap phan, qua chat so
                # voi do chinh xac 6 CHU SO CO NGHIA thuc te cua du lieu).
                tol = max(1e-5, 2e-5 * (abs(s) + abs(e)))
                self.assertLess(
                    abs(length - (e - s)), tol,
                    "{}: len_cm={} khong khop end-start={} (tol={})".format(label, length, e - s, tol),
                )
                g = int(row["generation"])
                sum_len_by_gen[g] = sum_len_by_gen.get(g, 0.0) + length

            total_chrom_len = self.REPLICATES * self.N_CHROM * self.CHROM_LENGTH_CM * 2.0
            for row in stats_rows:
                g = int(row["generation"])
                observed_sum = sum_len_by_gen.get(g, 0.0)
                expected_sum = float(row["mean_ancestry_fraction"]) * total_chrom_len
                tol = max(1e-2, expected_sum * 1e-3)  # dung sai lam tron 6 chu so + tich luy nhieu dong
                self.assertLess(
                    abs(observed_sum - expected_sum), tol,
                    "{} doi {}: tong len_cm={} khac mean_ancestry_fraction*tong do dai={} (tol={})".format(
                        label, g, observed_sum, expected_sum, tol
                    ),
                )


class TestSimulateReproducibilityAndGenerationZero(TempDirMixin, unittest.TestCase):
    def test_same_seed_gives_identical_output_across_two_runs(self):
        common_args = [
            "--generations", "4",
            "--mode", "isolated",
            "--replicates", "15",
            "--n-chrom", "3",
            "--chrom-length-cm", "80",
            "--seed", "123",
        ]
        out1_t, out1_s = self.path("run1_tracts.tsv"), self.path("run1_stats.tsv")
        out2_t, out2_s = self.path("run2_tracts.tsv"), self.path("run2_stats.tsv")
        rc1 = ds.main(["simulate"] + common_args + ["--out", out1_t, "--stats", out1_s])
        rc2 = ds.main(["simulate"] + common_args + ["--out", out2_t, "--stats", out2_s])
        self.assertEqual(rc1, 0)
        self.assertEqual(rc2, 0)
        self.assertEqual(read_text(out1_t), read_text(out2_t))
        self.assertEqual(read_text(out1_s), read_text(out2_s))

    def test_different_seed_generally_changes_output(self):
        base_args = [
            "--generations", "4",
            "--mode", "isolated",
            "--replicates", "15",
            "--n-chrom", "3",
            "--chrom-length-cm", "80",
        ]
        out1_t, out1_s = self.path("s1_tracts.tsv"), self.path("s1_stats.tsv")
        out2_t, out2_s = self.path("s2_tracts.tsv"), self.path("s2_stats.tsv")
        ds.main(["simulate"] + base_args + ["--seed", "1", "--out", out1_t, "--stats", out1_s])
        ds.main(["simulate"] + base_args + ["--seed", "2", "--out", out2_t, "--stats", out2_s])
        self.assertNotEqual(read_text(out1_t), read_text(out2_t))

    def test_generations_zero_returns_exactly_init_fraction(self):
        for init_fraction in (0.5, 0.3, 0.7):
            out_t, out_s = self.path("g0_tracts_{}.tsv".format(init_fraction)), self.path(
                "g0_stats_{}.tsv".format(init_fraction)
            )
            rc = ds.main(
                [
                    "simulate",
                    "--generations", "0",
                    "--mode", "backcross",
                    "--init-fraction", str(init_fraction),
                    "--replicates", "10",
                    "--n-chrom", "4",
                    "--chrom-length-cm", "100",
                    "--seed", "1",
                    "--out", out_t,
                    "--stats", out_s,
                ]
            )
            self.assertEqual(rc, 0)
            _, rows = tsv_rows(read_text(out_s))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["generation"], "0")
            self.assertAlmostEqual(float(rows[0]["mean_ancestry_fraction"]), init_fraction, places=6)
            self.assertEqual(float(rows[0]["sd_ancestry_fraction"]), 0.0)

    def test_invalid_arguments_raise_system_exit(self):
        base = ["--mode", "backcross", "--out", self.path("x.tsv"), "--stats", self.path("y.tsv")]
        with self.assertRaises(SystemExit):
            ds.main(["simulate", "--generations", "-1"] + base)
        with self.assertRaises(SystemExit):
            ds.main(["simulate", "--generations", "2", "--init-fraction", "1.5"] + base)
        with self.assertRaises(SystemExit):
            ds.main(["simulate", "--generations", "2", "--replicates", "0"] + base)
        with self.assertRaises(SystemExit):
            ds.main(["simulate", "--generations", "2", "--n-chrom", "0"] + base)
        with self.assertRaises(SystemExit):
            ds.main(["simulate", "--generations", "2", "--chrom-length-cm", "0"] + base)


# ---------------------------------------------------------------------------
# detect (qua CLI, main()) -- du lieu nhung (sim_stats.tsv nho, tu viet tay)
# ---------------------------------------------------------------------------


class TestDetectCommand(TempDirMixin, unittest.TestCase):
    def _write_stats(self, path, rows):
        lines = ["\t".join(ds.STATS_HEADER)]
        for r in rows:
            lines.append("\t".join(str(r[k]) for k in ds.STATS_HEADER))
        write_text(path, "\n".join(lines) + "\n")

    def test_thresholds_applied_with_inclusive_boundary(self):
        stats_path = self.path("stats.tsv")
        rows = [
            {  # doi 0: ca 2 dieu kien vua dung bang nguong -> ok (quy uoc >=)
                "generation": 0, "mean_ancestry_fraction": 0.5, "sd_ancestry_fraction": 0.0,
                "n_tracts": 10, "mean_tract_cm": 20.0, "median_tract_cm": 20.0, "max_tract_cm": 10.0,
            },
            {  # doi 1: fraction duoi nguong -> fraction_ok=0 -> detectable=0
                "generation": 1, "mean_ancestry_fraction": 0.009, "sd_ancestry_fraction": 0.01,
                "n_tracts": 5, "mean_tract_cm": 8.0, "median_tract_cm": 8.0, "max_tract_cm": 12.0,
            },
            {  # doi 2: tract duoi nguong -> tract_ok=0 -> detectable=0
                "generation": 2, "mean_ancestry_fraction": 0.5, "sd_ancestry_fraction": 0.01,
                "n_tracts": 5, "mean_tract_cm": 4.0, "median_tract_cm": 4.0, "max_tract_cm": 9.999,
            },
            {  # doi 3: ca 2 deu dat ro -> detectable=1
                "generation": 3, "mean_ancestry_fraction": 0.8, "sd_ancestry_fraction": 0.01,
                "n_tracts": 5, "mean_tract_cm": 15.0, "median_tract_cm": 15.0, "max_tract_cm": 20.0,
            },
        ]
        self._write_stats(stats_path, rows)
        out_path = self.path("detect.tsv")
        rc = ds.main(
            [
                "detect", "--stats", stats_path,
                "--min-fraction", "0.01", "--min-tract-cm", "10",
                "--out", out_path,
            ]
        )
        self.assertEqual(rc, 0)
        header, drows = tsv_rows(read_text(out_path))
        self.assertEqual(header, ds.DETECT_HEADER)
        by_gen = {int(r["generation"]): r for r in drows}

        self.assertEqual(by_gen[0]["fraction_ok"], "1")
        self.assertEqual(by_gen[0]["tract_ok"], "1")
        self.assertEqual(by_gen[0]["detectable"], "1")

        self.assertEqual(by_gen[1]["fraction_ok"], "0")
        self.assertEqual(by_gen[1]["detectable"], "0")

        self.assertEqual(by_gen[2]["tract_ok"], "0")
        self.assertEqual(by_gen[2]["detectable"], "0")

        self.assertEqual(by_gen[3]["detectable"], "1")

    def test_missing_columns_raises(self):
        bad_path = self.path("bad_stats.tsv")
        write_text(bad_path, "generation\tsome_other_col\n0\t1\n")
        with self.assertRaises(SystemExit):
            ds.main(
                [
                    "detect", "--stats", bad_path,
                    "--min-fraction", "0.01", "--min-tract-cm", "5",
                    "--out", self.path("out.tsv"),
                ]
            )


# ---------------------------------------------------------------------------
# report (qua CLI, main()) -- du lieu nhung, 2 kich ban tuong phan ro (giam
# manh kieu backcross vs phang kieu isolated) de kiem cau dien giai TU DONG
# doc dung xu huong, khong chi kiem co chay duoc.
# ---------------------------------------------------------------------------


class TestReportCommand(TempDirMixin, unittest.TestCase):
    def _write_stats(self, path, rows):
        lines = ["\t".join(ds.STATS_HEADER)]
        for r in rows:
            lines.append("\t".join(str(r[k]) for k in ds.STATS_HEADER))
        write_text(path, "\n".join(lines) + "\n")

    def _write_detect(self, path, rows):
        lines = ["\t".join(ds.DETECT_HEADER)]
        for r in rows:
            lines.append("\t".join(str(r[k]) for k in ds.DETECT_HEADER))
        write_text(path, "\n".join(lines) + "\n")

    def _stats_row(self, g, frac, mean_tract):
        return {
            "generation": g, "mean_ancestry_fraction": frac, "sd_ancestry_fraction": 0.01,
            "n_tracts": 10, "mean_tract_cm": mean_tract, "median_tract_cm": mean_tract,
            "max_tract_cm": mean_tract * 2,
        }

    def test_report_contains_mandatory_warning_and_mode_explanation(self):
        stats_path = self.path("stats.tsv")
        detect_path = self.path("detect.tsv")
        self._write_stats(stats_path, [self._stats_row(0, 0.5, 50.0), self._stats_row(1, 0.25, 30.0)])
        self._write_detect(
            detect_path,
            [
                {"generation": 0, "fraction_ok": 1, "tract_ok": 1, "detectable": 1},
                {"generation": 1, "fraction_ok": 1, "tract_ok": 1, "detectable": 1},
            ],
        )
        out_path = self.path("report.md")
        rc = ds.main(["report", "--stats", stats_path, "--detect", detect_path, "--out", out_path])
        self.assertEqual(rc, 0)
        report = read_text(out_path)
        flat_report = normalize_ws(report)  # Markdown gap xuong dong mem thanh dau cach khi render
        self.assertIn("CANH BAO BAT BUOC", report)
        self.assertIn("Chon loc", report)
        self.assertIn("Thay doi kich thuoc quan the", report)
        self.assertIn("Bien thien toc do tai to hop", report)
        self.assertIn("Ban do di truyen dong nhat", report)
        self.assertIn("BAC DO LON", report)
        self.assertIn("backcross", report)
        self.assertIn("isolated", report)
        self.assertIn("THAM SO NGUOI DUNG DAT", flat_report)
        self.assertIn("| generation | mean_ancestry_fraction |", report)

    def test_report_detects_backcross_like_decay_pattern(self):
        stats_path = self.path("stats.tsv")
        detect_path = self.path("detect.tsv")
        rows = [self._stats_row(g, 0.5 * (0.5 ** g), 50.0 / (g + 1)) for g in range(6)]
        self._write_stats(stats_path, rows)
        self._write_detect(
            detect_path,
            [{"generation": g, "fraction_ok": (1 if g < 3 else 0), "tract_ok": 1, "detectable": (1 if g < 3 else 0)} for g in range(6)],
        )
        out_path = self.path("report.md")
        ds.main(["report", "--stats", stats_path, "--detect", detect_path, "--out", out_path])
        report = read_text(out_path)
        self.assertIn("giam manh, don dieu", report)
        self.assertIn("phu hop kieu 'lai nguoc' (backcross)", report)
        self.assertIn("con phat hien duoc toi da den doi 2", report)

    def test_report_detects_isolated_like_flat_pattern(self):
        stats_path = self.path("stats.tsv")
        detect_path = self.path("detect.tsv")
        # dao dong nho quanh 0.5, khong giam he thong
        fracs = [0.5, 0.52, 0.49, 0.51, 0.5, 0.505]
        rows = [self._stats_row(g, fracs[g], 40.0 - g) for g in range(6)]
        self._write_stats(stats_path, rows)
        self._write_detect(
            detect_path,
            [{"generation": g, "fraction_ok": 1, "tract_ok": 1, "detectable": 1} for g in range(6)],
        )
        out_path = self.path("report.md")
        ds.main(["report", "--stats", stats_path, "--detect", detect_path, "--out", out_path])
        report = read_text(out_path)
        self.assertIn("khong giam he thong", report)
        self.assertIn("giao phoi noi bo", report)

    def test_report_empty_stats_raises(self):
        empty_path = self.path("empty.tsv")
        write_text(empty_path, "")
        with self.assertRaises(SystemExit):
            ds.main(
                [
                    "report", "--stats", empty_path, "--detect", empty_path,
                    "--out", self.path("report.md"),
                ]
            )


# ---------------------------------------------------------------------------
# --help (top-level va tung lenh con)
# ---------------------------------------------------------------------------


class TestHelp(unittest.TestCase):
    def _assert_help_exits_zero(self, argv):
        with self.assertRaises(SystemExit) as ctx:
            ds.main(argv)
        self.assertEqual(ctx.exception.code, 0)

    def test_top_level_help(self):
        self._assert_help_exits_zero(["--help"])

    def test_simulate_help(self):
        self._assert_help_exits_zero(["simulate", "--help"])

    def test_detect_help(self):
        self._assert_help_exits_zero(["detect", "--help"])

    def test_report_help(self):
        self._assert_help_exits_zero(["report", "--help"])


if __name__ == "__main__":
    unittest.main()
