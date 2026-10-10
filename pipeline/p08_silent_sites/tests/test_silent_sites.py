"""Unit test cho silent_sites.py -- unittest thuan, du lieu nhung, khong mang.

Chay: python -B -m unittest discover -s pipeline/p08_silent_sites/tests -v
"""

import gzip
import io
import json
import math
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import silent_sites as ss  # noqa: E402

# Back-translation uu tien ho codon 4-fold de co du vi tri cam
_BACK = {
    "L": "CTG", "V": "GTT", "S": "TCT", "P": "CCT", "T": "ACT", "A": "GCT", "R": "CGT", "G": "GGT",
    "M": "ATG", "K": "AAA", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAG", "E": "GAA", "H": "CAT",
    "I": "ATT", "F": "TTT", "W": "TGG", "Y": "TAT",
}
_ROT = {"T": "C", "C": "A", "A": "G", "G": "T"}


def make_protein(n=120, seed=3):
    rnd = random.Random(seed)
    return "M" + "".join(rnd.choice(list(_BACK)) for _ in range(n - 1))


def backtranslate(prot):
    return "".join(_BACK[a] for a in prot)


def rotate_fourfold(cds):
    """Doi base 3 cua MOI codon 4-fold (dong nghia) -> protein giu nguyen."""
    out = []
    for i in range(0, len(cds), 3):
        c = cds[i : i + 3]
        if c[:2] in ss.FOURFOLD_PREFIXES:
            c = c[:2] + _ROT[c[2]]
        out.append(c)
    return "".join(out)


def randomize_fourfold(cds, seed=11):
    rnd = random.Random(seed)
    out = []
    for i in range(0, len(cds), 3):
        c = cds[i : i + 3]
        if c[:2] in ss.FOURFOLD_PREFIXES:
            c = c[:2] + rnd.choice("ACGT")
        out.append(c)
    return "".join(out)


def read_text(path):
    with open(path, "r", newline="", encoding="utf-8") as f:
        return f.read()


def write_text(path, content):
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(content)


def tsv_rows(text):
    lines = [ln for ln in text.splitlines() if ln != ""]
    header = lines[0].split("\t")
    return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


class TempDirMixin:
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ss_test_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)


# ---------------------------------------------------------------------------


class TestCodons(unittest.TestCase):
    def test_translate(self):
        self.assertEqual(ss.translate("ATGTTTTAA"), "MF*")
        self.assertEqual(ss.translate("ATGNNN"), "MX")

    def test_fourfold_pair(self):
        self.assertTrue(ss.is_fourfold_pair("CTT", "CTG"))  # Leu, cung prefix 4-fold
        self.assertTrue(ss.is_fourfold_pair("GGA", "GGA"))
        self.assertFalse(ss.is_fourfold_pair("CTT", "GTT"))  # prefix khac
        self.assertFalse(ss.is_fourfold_pair("ATT", "ATC"))  # Ile: khong phai 4-fold
        self.assertFalse(ss.is_fourfold_pair("TTA", "TTG"))  # Leu TTR: 2-fold
        self.assertFalse(ss.is_fourfold_pair("CTN", "CTG"))  # ky tu mo ho

    def test_jukes_cantor(self):
        self.assertEqual(ss.jukes_cantor(0.0), 0.0)
        self.assertTrue(math.isinf(ss.jukes_cantor(0.75)))
        self.assertTrue(math.isinf(ss.jukes_cantor(0.9)))
        self.assertAlmostEqual(ss.jukes_cantor(0.3), -0.75 * math.log(0.6), places=9)

    def test_strip_stop_and_check(self):
        self.assertEqual(ss.strip_stop("ATGTTTTAA"), "ATGTTT")
        self.assertEqual(ss.strip_stop("ATGTTT"), "ATGTTT")
        self.assertEqual(ss.check_cds("ATGTTTT"), (False, "len_not_multiple_of_3"))
        self.assertEqual(ss.check_cds("ATGTTT"), (False, "too_short"))
        good = backtranslate(make_protein(40))
        self.assertEqual(ss.check_cds(good), (True, ""))
        bad = good[:30] + "TAA" + good[33:]
        self.assertEqual(ss.check_cds(bad), (False, "internal_stop"))

    def test_blosum62_symmetric_positive_diagonal(self):
        for a in ss._B62_ORDER:
            self.assertGreater(ss.BLOSUM62[(a, a)], 0)
            for b in ss._B62_ORDER:
                self.assertEqual(ss.BLOSUM62[(a, b)], ss.BLOSUM62[(b, a)])
        self.assertEqual(ss.b62("W", "W"), 11)
        self.assertEqual(ss.b62("X", "A"), -1)


class TestAlign(unittest.TestCase):
    def test_identical_no_gaps(self):
        p = make_protein(60)
        a1, a2, score = ss.align_global(p, p)
        self.assertEqual(a1, p)
        self.assertEqual(a2, p)
        self.assertEqual(score, sum(ss.b62(x, x) for x in p))

    def test_insertion_placed_as_gap(self):
        a1, a2, _ = ss.align_global("MKTAYIAKQRQISFVK", "MKTAYQRQISFVK")
        self.assertEqual(a1, "MKTAYIAKQRQISFVK")
        self.assertEqual(a2.count("-"), 3)
        self.assertEqual(a2.replace("-", ""), "MKTAYQRQISFVK")
        b1, b2, _ = ss.align_global("MKTAYQRQISFVK", "MKTAYIAKQRQISFVK")
        self.assertEqual(b1.count("-"), 3)

    def test_unknown_residue_does_not_crash(self):
        a1, a2, _ = ss.align_global("MKXAY", "MKTAY")
        self.assertEqual(len(a1), len(a2))


class TestSiteStats(unittest.TestCase):
    def setUp(self):
        self.prot = make_protein(120)
        self.cds = backtranslate(self.prot)

    def _stats(self, cds2):
        p2 = ss.translate(cds2)
        a1, a2, _ = ss.align_global(self.prot, p2)
        cols = ss.codon_columns(a1, a2, self.cds, cds2)
        return cols, ss.site_stats(cols)

    def test_identical(self):
        cols, s = self._stats(self.cds)
        self.assertEqual(s["n_codon"], 120)
        self.assertEqual(s["aa_ident"], 1.0)
        self.assertEqual(s["ident_p3"], 1.0)
        self.assertGreaterEqual(s["n_4fold"], 20)
        self.assertEqual(s["ident_4fold"], 1.0)
        self.assertEqual(s["jc_4fold"], 0.0)

    def test_all_fourfold_changed_keeps_protein(self):
        cds2 = rotate_fourfold(self.cds)
        self.assertEqual(ss.translate(cds2), self.prot)
        cols, s = self._stats(cds2)
        self.assertEqual(s["aa_ident"], 1.0)
        self.assertEqual(s["ident_p1"], 1.0)
        self.assertEqual(s["ident_p2"], 1.0)
        self.assertEqual(s["ident_4fold"], 0.0)
        self.assertTrue(math.isinf(s["jc_4fold"]))
        self.assertLess(s["ident_p3"], 1.0)

    def test_randomized_fourfold_near_chance(self):
        cds2 = randomize_fourfold(self.cds)
        self.assertEqual(ss.translate(cds2), self.prot)
        cols, s = self._stats(cds2)
        self.assertEqual(s["aa_ident"], 1.0)
        self.assertLess(s["ident_4fold"], 0.6)
        self.assertGreater(s["chance_4fold"], 0.15)
        self.assertLess(s["chance_4fold"], 0.6)
        self.assertGreater(s["jc_4fold"], 0.3)

    def test_core_mask_drops_divergent_block(self):
        # thay 30 aa giua (vi tri 45-74) bang axit amin TUONG TU (BLOSUM62 >= -1)
        # de can khong mo gap nhung khong con giong het -> vung do khong la loi.
        # Cua so +/-10: cot 44 con 11/21 giong (loi), cot 45 chi 10/21 (ngoai loi).
        similar = {"A": "S", "S": "T", "T": "S", "L": "I", "I": "V", "V": "I", "K": "R", "R": "K",
                   "D": "E", "E": "D", "N": "D", "Q": "E", "F": "Y", "Y": "F", "G": "A", "P": "A",
                   "M": "L", "W": "F", "C": "S", "H": "Y"}
        p2 = self.prot[:45] + "".join(similar[c] for c in self.prot[45:75]) + self.prot[75:]
        cds2 = backtranslate(p2)
        a1, a2, _ = ss.align_global(self.prot, p2)
        self.assertNotIn("-", a1 + a2)
        cols = ss.codon_columns(a1, a2, self.cds, cds2)
        self.assertEqual(len(cols), 120)
        mask = ss.core_mask(cols)
        self.assertTrue(all(mask[:45]))
        self.assertFalse(any(mask[45:75]))
        self.assertTrue(all(mask[75:]))
        s_core = ss.site_stats(cols, mask)
        s_all = ss.site_stats(cols)
        self.assertGreater(s_core["aa_ident"], s_all["aa_ident"])
        self.assertLess(s_core["n_codon"], s_all["n_codon"])

    def test_codon_columns_skip_gaps(self):
        cols = ss.codon_columns("MK-AY", "MKTAY", "ATGAAAGCTTAT", "ATGAAAACTGCTTAT")
        self.assertEqual([c[0] for c in cols], ["ATG", "AAA", "GCT", "TAT"])
        self.assertEqual([c[1] for c in cols], ["ATG", "AAA", "GCT", "TAT"])

    def test_gene_class(self):
        self.assertEqual(ss.gene_class("HOXA1"), "hox")
        self.assertEqual(ss.gene_class("hoxd13"), "hox")
        self.assertEqual(ss.gene_class("PAX6"), "toolkit")
        self.assertEqual(ss.gene_class("ACTB"), "control")
        self.assertEqual(ss.gene_class("ACTB", {"ACTB": "toolkit"}), "toolkit")


# ---------------------------------------------------------------------------


def _ncbi_header(acc, pid, gene):
    return "lcl|{a}_cds_{p}_1 [gene={g}] [protein_id={p}] [protein=test]".format(a=acc, p=pid, g=gene)


class TestCompareCLI(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.prot = make_protein(120, seed=5)
        self.cds = backtranslate(self.prot)
        self.human = self.path("human.fna")
        with open(self.human, "w", newline="", encoding="utf-8") as f:
            ss.write_fasta_record(f, "HOXA1|human|NM_1.1|NP_1.1 [gene=HOXA1] [protein_id=NP_1.1]", self.cds + "TAA")
            ss.write_fasta_record(f, "ACTB|human|NM_2.1|NP_2.1 [gene=ACTB] [protein_id=NP_2.1]", self.cds + "TGA")
            ss.write_fasta_record(f, "PAX6|human|NM_3.1|NP_3.1 [gene=PAX6] [protein_id=NP_3.1]", self.cds)
            ss.write_fasta_record(f, "SOX9|human|NM_4.1|NP_4.1 [gene=SOX9] [protein_id=NP_4.1]", self.cds)
        self.bulbul = self.path("bulbul_cds.fna.gz")
        with gzip.open(self.bulbul, "wt", encoding="utf-8", newline="") as f:
            # Hoxa1: 2 isoform, ban dai hon la ban vi tri cam ngau nhien
            ss.write_fasta_record(f, _ncbi_header("VW1.1", "XP_11.1", "Hoxa1"), randomize_fourfold(self.cds) + "TAA")
            ss.write_fasta_record(f, _ncbi_header("VW1.1", "XP_12.1", "Hoxa1"), self.cds[:90] + "TAA")
            # Actb: ban sao giong het (kich ban "chep gan day")
            ss.write_fasta_record(f, _ncbi_header("VW2.1", "XP_21.1", "Actb"), self.cds + "TAA")
            # Pax6: codon dung noi -> bo
            ss.write_fasta_record(f, _ncbi_header("VW3.1", "XP_31.1", "Pax6"), self.cds[:60] + "TAG" + self.cds[63:])
            # gene khong lien quan
            ss.write_fasta_record(f, _ncbi_header("VW4.1", "XP_41.1", "Efcc1"), self.cds)

    def test_compare_end_to_end(self):
        out = self.path("stats.tsv")
        aln = self.path("aln.txt")
        rc = ss.main(["compare", "--human", self.human, "--bulbul", self.bulbul, "--out", out, "--alignments", aln])
        self.assertEqual(rc, 0)
        h, rows = tsv_rows(read_text(out))
        self.assertEqual(h, ss.STATS_HEADER)
        by = {r["gene"]: r for r in rows}
        self.assertEqual(set(by), {"HOXA1", "ACTB", "PAX6", "SOX9"})
        self.assertEqual([r["gene"] for r in rows], ["HOXA1", "PAX6", "SOX9", "ACTB"])  # hox, toolkit, control
        hox = by["HOXA1"]
        self.assertEqual(hox["class"], "hox")
        self.assertEqual(hox["bulbul_id"], "lcl|VW1.1_cds_XP_11.1_1")  # isoform dai nhat
        self.assertEqual(hox["status"], "ok")
        self.assertEqual(float(hox["aa_ident_core"]), 1.0)
        self.assertLess(float(hox["ident_4fold_core"]), 0.6)
        actb = by["ACTB"]
        self.assertEqual(actb["class"], "control")
        self.assertEqual(float(actb["ident_4fold_core"]), 1.0)
        self.assertEqual(float(actb["jc_4fold_core"]), 0.0)
        self.assertEqual(by["PAX6"]["status"], "bulbul_internal_stop")
        self.assertEqual(by["SOX9"]["status"], "no_bulbul_cds")
        txt = read_text(aln)
        self.assertIn("## HOXA1", txt)
        self.assertIn("core   ", txt)
        self.assertNotIn("## PAX6", txt)

    def test_report_from_compare(self):
        out = self.path("stats.tsv")
        ss.main(["compare", "--human", self.human, "--bulbul", self.bulbul, "--out", out])
        rep = self.path("r.md")
        self.assertEqual(ss.main(["report", "--tsv", out, "--out", rep, "--date", "2026-10-10"]), 0)
        txt = read_text(rep)
        self.assertIn("## 1. Gene Hox", txt)
        self.assertIn("| HOXA1 |", txt)
        self.assertIn("khớp thừa kế lâu đời", txt)  # HOXA1 ngau nhien
        self.assertIn("dấu hiệu chép/ghép gần đây", txt)  # ACTB giong het
        self.assertIn("không có CDS chào mào mang tên gene này", txt)  # SOX9
        self.assertIn("không xét (bulbul_internal_stop)", txt)  # PAX6
        self.assertIn("Cách đọc — BẮT BUỘC", txt)
        self.assertIn("Cấm viết \"đã chứng minh\"", txt)
        self.assertIn("## 4. So Hox với đối chứng", txt)

    def test_verdict_rules(self):
        base = {"status": "ok", "ident_4fold_core": "0.95", "chance_4fold_core": "0.3"}
        self.assertIn("chép/ghép", ss.verdict(base))
        self.assertIn("ngẫu nhiên →", ss.verdict({"status": "ok", "ident_4fold_core": "0.35", "chance_4fold_core": "0.3"}))
        self.assertIn("ràng buộc codon", ss.verdict({"status": "ok", "ident_4fold_core": "0.6", "chance_4fold_core": "0.3"}))
        self.assertIn("dưới ngẫu nhiên", ss.verdict({"status": "ok", "ident_4fold_core": "0.05", "chance_4fold_core": "0.3"}))
        self.assertIn("không xét", ss.verdict({"status": "few_4fold_sites", "ident_4fold_core": "1", "chance_4fold_core": "0.3"}))


# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, body):
        self._b = body.encode("utf-8")

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _qs(url):
    import urllib.parse as up

    return {k: v[0] for k, v in up.parse_qs(up.urlparse(url).query).items()}


REGISTRY_HEADER = "dataset_id,name,source_url,version_or_release,downloaded_on,local_path,size_bytes,sha256,license,notes\n"


class TestFetchHuman(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.calls = []
        self._iv = ss.REQUEST_INTERVAL
        ss.REQUEST_INTERVAL = 0.0
        self.addCleanup(setattr, ss, "REQUEST_INTERVAL", self._iv)
        prot = make_protein(60, seed=9)
        self.cds = backtranslate(prot)

        def fake(req, timeout=60):
            url = req.full_url
            self.calls.append(url)
            q = _qs(url)
            if "esearch.fcgi" in url:
                self.assertEqual(q["db"], "nuccore")
                term = q["term"]
                ids = ["1", "2"] if term.startswith("HOXA1") else []
                return _FakeResp(json.dumps({"esearchresult": {"idlist": ids}}))
            if "efetch.fcgi" in url:
                self.assertEqual(q["rettype"], "fasta_cds_na")
                fa = ">lcl|NM_100.1_cds_NP_100.1_1 [gene=HOXA1] [protein_id=NP_100.1] [protein=homeobox A1 isoform 1]\n" + self.cds + "TAA\n"
                fa += ">lcl|NM_101.2_cds_NP_101.1_1 [gene=HOXA1] [protein_id=NP_101.1] [protein=homeobox A1 isoform 2]\n" + self.cds[:120] + "TAA\n"
                fa += ">lcl|NM_200.1_cds_NP_200.1_1 [gene=OTHER] [protein_id=NP_200.1] [protein=wrong gene]\n" + self.cds + self.cds + "TAA\n"
                return _FakeResp(fa)
            raise AssertionError(url)

        self._orig = ss._URLOPEN
        ss._URLOPEN = fake
        self.addCleanup(setattr, ss, "_URLOPEN", self._orig)

    def test_fetch_picks_longest_right_gene(self):
        out = self.path("h.fna")
        man = self.path("m.tsv")
        reg = self.path("registry.csv")
        write_text(reg, REGISTRY_HEADER + "D01-x,old,u,v,2026-09-07,p,1,abc,pub,n\n")
        rc = ss.main(["fetch-human", "--genes", "HOXA1,HOXB13", "--email", "t@example.org", "--out", out, "--manifest", man, "--registry", reg, "--date", "2026-10-10"])
        self.assertEqual(rc, 0)
        recs = list(ss.iter_fasta(io.StringIO(read_text(out))))
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0][0], "HOXA1|human|NM_100.1|NP_100.1")
        self.assertEqual(recs[0][2], self.cds + "TAA")
        _, rows = tsv_rows(read_text(man))
        by = {r["gene"]: r for r in rows}
        self.assertEqual(by["HOXA1"]["status"], "ok")
        self.assertEqual(by["HOXA1"]["n_candidates"], "2")
        self.assertEqual(by["HOXA1"]["class"], "hox")
        self.assertEqual(by["HOXB13"]["status"], "not_found")
        self.assertIn("biomol_mrna[PROP]", by["HOXA1"]["search_term"])
        lines = read_text(reg).splitlines()
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[2].startswith("D10-b11-human-cds-"))
        self.assertIn(ss.sha256_file(out), lines[2])
        # idempotent
        ss.main(["fetch-human", "--genes", "HOXA1", "--email", "t@example.org", "--out", out, "--manifest", man, "--registry", reg, "--date", "2026-10-10"])
        self.assertEqual(len(read_text(reg).splitlines()), 3)

    def test_requires_email(self):
        with self.assertRaises(SystemExit):
            ss.main(["fetch-human", "--email", "x", "--out", self.path("a"), "--manifest", self.path("b")])


class TestHelp(unittest.TestCase):
    def test_help_exit_zero(self):
        for argv in (["--help"], ["fetch-human", "--help"], ["compare", "--help"], ["report", "--help"]):
            with self.assertRaises(SystemExit) as cm:
                ss.main(argv)
            self.assertEqual(cm.exception.code, 0, argv)


if __name__ == "__main__":
    unittest.main()
