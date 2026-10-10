"""Unit test cho tooth_scan.py -- unittest thuan, khong pytest, du lieu nhung
ngay trong file (khong doc/ghi ngoai thu muc tam cua test, khong mang).

Chay: python -B -m unittest discover -s pipeline/p07_tooth_genes/tests -v

Bang HSP nhung duoi day (HSP_ROWS) la DAU RA THAT cua tblastn 2.12.0+ chay tren
bo gene tong hop sinh boi make_synthetic_genome() (seed co dinh), gom 4
scaffold: nguyen ven / co codon dung + lech khung (mach +) / cung locus do
dao nguoc (mach -) / 2 exon cach 1 intron 400 nt. Test tich hop cuoi file
chay lai makeblastdb + tblastn THAT neu 2 lenh do co tren PATH (bo qua neu
khong), de bao dam chuoi -outfmt OUTFMT van doc duoc o phien ban BLAST khac.
"""

import gzip
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tooth_scan as ts  # noqa: E402


# ---------------------------------------------------------------------------
# Du lieu tong hop (cung bo sinh ra HSP_ROWS)
# ---------------------------------------------------------------------------

_BACK = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAG", "E": "GAA",
    "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT",
    "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
}


def make_synthetic_genome(seed=7):
    """-> (protein 180 aa, {scaffold: seq}). Loci:
    scaf_intact   : CDS nguyen ven (+)
    scaf_lesion   : codon 91 -> TAA (codon dung) + xoa 1 nt o codon 141 (lech khung), mach +
    scaf_lesion_rc: dao nguoc bo sung scaf_lesion (mach -)
    scaf_intron   : exon1 (100 aa) + intron 400 nt (GT..AG) + exon2 (80 aa)"""
    rnd = random.Random(seed)
    aas = list(_BACK)
    prot = "M" + "".join(rnd.choice(aas) for _ in range(179))
    cds = "".join(_BACK[a] for a in prot)

    def rand_nt(n):
        return "".join(rnd.choice("ACGT") for _ in range(n))

    cds_stop = cds[:270] + "TAA" + cds[273:]
    cds_fs = cds_stop[:420] + cds_stop[421:]
    ex1, ex2 = cds[:300], cds[300:]
    intron = "GT" + rand_nt(396) + "AG"
    genome = {
        "scaf_intact": rand_nt(500) + cds + rand_nt(500),
        "scaf_lesion": rand_nt(300) + cds_fs + rand_nt(300),
        "scaf_intron": rand_nt(200) + ex1 + intron + ex2 + rand_nt(200),
    }
    genome["scaf_lesion_rc"] = ts.revcomp(genome["scaf_lesion"])
    return prot, genome


HSP_ROWS = """ENAM_test	scaf_intact	100.000	180	0	0	1	180	501	1040	2.13e-135	378	180	1540	3	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL
ENAM_test	scaf_lesion_rc	99.286	140	1	0	1	140	839	420	4.19e-103	291	180	1139	-1	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGC	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNW*LLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGC
ENAM_test	scaf_lesion_rc	87.500	48	6	0	133	180	444	301	5.80e-25	88.6	180	1139	-3	FWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL	FLAF*NGLYQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL
ENAM_test	scaf_lesion	99.286	140	1	0	1	140	301	720	4.19e-103	291	180	1139	1	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGC	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNW*LLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGC
ENAM_test	scaf_lesion	87.500	48	6	0	133	180	696	839	5.80e-25	88.6	180	1139	3	FWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL	FLAF*NGLYQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL
ENAM_test	scaf_intron	86.066	122	16	1	1	122	201	563	2.50e-72	214	180	1340	3	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICG	MLCMRNWDKYRTERNFFNGNWFRYDGYRYYMRGRWCIFCWDYIWQDYYEKDWNYRVESWFLPYPKIGQGNYITSLPIVNDTFQLCSFRNWYLLKVSYPNN-VNRTRVILTGHAEARPPEVRG
ENAM_test	scaf_intron	77.982	109	19	2	72	180	829	1140	9.10e-57	174	180	1340	1	ITSLPIVNDTFQLCSFRNWYLLKVSYPNNHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL	ITND*LINEPFMTRAY-DWFTIVC----QHSNRIYPIMKAPKQVDSREICGMMSNQPMWHCFWHFKMGCNQCGGASYQHIACFWKVYLCTVRPWMMMMDSMRENEPQDL
"""


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
        self.tmp = tempfile.mkdtemp(prefix="ts_test_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)

    def hsp_file(self):
        p = self.path("hsps.tsv")
        write_text(p, HSP_ROWS)
        return p


# ---------------------------------------------------------------------------
# Dich ma
# ---------------------------------------------------------------------------


class TestTranslate(unittest.TestCase):
    def test_known_codons(self):
        self.assertEqual(ts.translate("ATGTTTTAA"), "MF*")
        self.assertEqual(ts.translate("ATGTAGTGA"), "M**")
        self.assertEqual(ts.translate("ATGNNNTTT"), "MXF")

    def test_frames(self):
        nt = "GATGTTTAAC"  # khung +2 (frame=1): ATG TTT AAC -> MFN
        self.assertEqual(ts.translate(nt, 1), "MFN")
        self.assertEqual(ts.translate(nt, 0), "DV*")

    def test_revcomp(self):
        self.assertEqual(ts.revcomp("ATGCN"), "NGCAT")
        self.assertEqual(ts.revcomp(ts.revcomp("ACGTTGCA")), "ACGTTGCA")

    def test_three_frame_report_marks_stops(self):
        txt = ts.three_frame_report("x", "ATGTAATTT")
        self.assertIn("khung +1: 1 codon dung tai nt 4", txt)
        self.assertIn("M*F", txt)


# ---------------------------------------------------------------------------
# read_hsps / chain_loci
# ---------------------------------------------------------------------------


class TestReadHsps(TempDirMixin, unittest.TestCase):
    def test_reads_all_rows_and_strand(self):
        hsps = ts.read_hsps(self.hsp_file())
        self.assertEqual(len(hsps), 7)
        by_scaf = {}
        for h in hsps:
            by_scaf.setdefault(h["sseqid"], []).append(h)
        self.assertEqual(set(by_scaf), {"scaf_intact", "scaf_lesion", "scaf_lesion_rc", "scaf_intron"})
        for h in by_scaf["scaf_lesion_rc"]:
            self.assertEqual(h["strand"], "-")
            self.assertLess(h["send"], h["sstart"])
            self.assertEqual(h["s_lo"], h["send"])
        for h in by_scaf["scaf_lesion"]:
            self.assertEqual(h["strand"], "+")
        self.assertIsInstance(hsps[0]["evalue"], float)
        self.assertIsInstance(hsps[0]["qstart"], int)

    def test_evalue_filter(self):
        all_h = ts.read_hsps(self.hsp_file())
        strict = ts.read_hsps(self.hsp_file(), max_evalue=1e-60)
        self.assertLess(len(strict), len(all_h))
        self.assertTrue(all(h["evalue"] <= 1e-60 for h in strict))

    def test_wrong_column_count_exits(self):
        p = self.path("bad.tsv")
        write_text(p, "a\tb\tc\n")
        with self.assertRaises(SystemExit) as cm:
            ts.read_hsps(p)
        self.assertIn("-outfmt", str(cm.exception))

    def test_comment_and_blank_lines_skipped(self):
        p = self.path("c.tsv")
        write_text(p, "# comment\n\n" + HSP_ROWS)
        self.assertEqual(len(ts.read_hsps(p)), 7)


class TestChainLoci(TempDirMixin, unittest.TestCase):
    def test_default_gap_gives_four_loci(self):
        loci = ts.chain_loci(ts.read_hsps(self.hsp_file()), max_locus_gap=100000)
        keys = sorted((l["sseqid"], l["strand"], len(l["hsps"])) for l in loci)
        self.assertEqual(
            keys,
            [("scaf_intact", "+", 1), ("scaf_intron", "+", 2), ("scaf_lesion", "+", 2), ("scaf_lesion_rc", "-", 2)],
        )

    def test_small_gap_splits_intron_locus(self):
        loci = ts.chain_loci(ts.read_hsps(self.hsp_file()), max_locus_gap=100)
        intron = [l for l in loci if l["sseqid"] == "scaf_intron"]
        self.assertEqual(len(intron), 2)  # 2 exon cach ~266 nt > 100 -> tach
        lesion = [l for l in loci if l["sseqid"] == "scaf_lesion"]
        self.assertEqual(len(lesion), 1)  # 2 HSP chong lan -> van 1 locus


# ---------------------------------------------------------------------------
# analyze_locus: codon dung, lech khung, ranh exon
# ---------------------------------------------------------------------------


def _locus(hsps, scaffold):
    loci = ts.chain_loci(hsps, 100000)
    return [l for l in loci if l["sseqid"] == scaffold][0]


class TestAnalyzeLocus(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.hsps = ts.read_hsps(self.hsp_file())

    def test_intact_locus_clean(self):
        s, les = ts.analyze_locus(_locus(self.hsps, "scaf_intact"))
        self.assertAlmostEqual(s["qcov"], 1.0)
        self.assertEqual(s["n_hsp"], 1)
        self.assertEqual((s["n_stop_raw"], s["n_stop_supported"], s["n_frameshift"], s["n_intron_like"]), (0, 0, 0, 0))
        self.assertEqual(les, [])
        self.assertEqual((s["s_start"], s["s_end"]), (501, 1040))

    def test_lesion_plus_strand(self):
        s, les = ts.analyze_locus(_locus(self.hsps, "scaf_lesion"))
        # codon 91 = TAA -> dung 1 '*' co cho dua, dung vi tri truy van 91,
        # toa do nt = 301 + 90*3 = 571 (CDS bat dau 301 1-based).
        self.assertEqual(s["n_stop_supported"], 1)
        stops = [l for l in les if l["type"] == "stop_codon" and l["supported"] == 1]
        self.assertEqual(len(stops), 1)
        self.assertEqual(stops[0]["qpos"], 91)
        self.assertEqual(stops[0]["s_nt"], 571)
        # '*' gia trong HSP yeu mo rong nguoc vao khung sai: co trong n_stop_raw,
        # KHONG duoc tinh (owned=0)
        self.assertEqual(s["n_stop_raw"], 2)
        spurious = [l for l in les if l["type"] == "stop_codon" and l["supported"] == 0]
        self.assertEqual(len(spurious), 1)
        self.assertEqual(spurious[0]["owned"], 0)
        # lech khung: khung 1 -> 3, subject gap am (chong lan), 1 ung vien
        self.assertEqual(s["n_frameshift"], 1)
        fs = [l for l in les if l["type"] == "frameshift"][0]
        self.assertEqual((fs["frame_from"], fs["frame_to"]), (1, 3))
        self.assertLessEqual(fs["gap_nt"], 40)
        self.assertEqual(s["n_intron_like"], 0)
        self.assertGreater(s["qcov"], 0.95)

    def test_lesion_minus_strand_mirrors_plus(self):
        s, les = ts.analyze_locus(_locus(self.hsps, "scaf_lesion_rc"))
        self.assertEqual(s["strand"], "-")
        self.assertEqual(s["n_stop_supported"], 1)
        self.assertEqual(s["n_frameshift"], 1)
        stop = [l for l in les if l["type"] == "stop_codon" and l["supported"] == 1][0]
        self.assertEqual(stop["qpos"], 91)
        # scaffold dai 1139; codon 571-573 tren mach + goc -> 567-569 tren
        # toa do thuan cua scaffold dao nguoc; s_nt = base thap nhat = 567
        self.assertEqual(stop["s_nt"], 567)
        fs = [l for l in les if l["type"] == "frameshift"][0]
        self.assertEqual((fs["frame_from"], fs["frame_to"]), (-1, -3))

    def test_intron_boundary_not_a_lesion(self):
        s, les = ts.analyze_locus(_locus(self.hsps, "scaf_intron"))
        self.assertEqual(s["n_intron_like"], 1)
        self.assertEqual(s["n_frameshift"], 0)
        # HSP 2 mo rong nguoc vao intron sinh 1 '*' gia -> raw 1, supported 0
        self.assertEqual(s["n_stop_raw"], 1)
        self.assertEqual(s["n_stop_supported"], 0)
        kinds = sorted(l["type"] for l in les)
        self.assertEqual(kinds, ["intron_like", "stop_codon"])
        self.assertGreater(s["qcov"], 0.95)

    def test_flank_support_requires_both_sides(self):
        # cot: (col, qpos, s_nt, qch, sch); '*' o giua, trai khop 3, phai khop 3
        cols = []
        for i in range(3):
            cols.append((i, i + 1, None, "A", "A"))
        cols.append((3, 4, None, "K", "*"))
        for i in range(4, 7):
            cols.append((i, i + 1, None, "A", "A"))
        self.assertTrue(ts._flank_supported(cols, 3, flank=3, min_ident=1.0))
        self.assertFalse(ts._flank_supported(cols, 3, flank=4, min_ident=1.0))
        cols2 = list(cols)
        cols2[5] = (5, 6, None, "A", "G")  # 1 mismatch ben phai -> 2/3
        self.assertFalse(ts._flank_supported(cols2, 3, flank=3, min_ident=0.9))
        self.assertTrue(ts._flank_supported(cols2, 3, flank=3, min_ident=0.6))

    def test_union_len(self):
        self.assertEqual(ts._union_len([(1, 10), (5, 12), (20, 25)]), 18)
        self.assertEqual(ts._union_len([]), 0)
        self.assertEqual(ts._union_len([(1, 3), (4, 6)]), 6)  # ke nhau -> noi


# ---------------------------------------------------------------------------
# CLI parse / extract / report
# ---------------------------------------------------------------------------

MANIFEST = (
    "qseqid\tgene\trole\taccession\torganism\tlength_aa\tn_candidates\tstatus\ttitle\tsearch_term\n"
    "ENAM_test\tENAM\ttarget\tXP_TEST.1\tAlligator test\t180\t1\tok\ttest protein\tENAM[Gene Name]\n"
    "AMBN_x\tAMBN\ttarget\t\tAlligator test\t0\t0\tnot_found\t\tAMBN[Gene Name]\n"
    "MEPE_x\tMEPE\tcontrol\tXP_CTRL.1\tAlligator test\t100\t1\tok\tctrl\tMEPE[Gene Name]\n"
)


class TestParseCLI(TempDirMixin, unittest.TestCase):
    def run_parse(self, extra=()):
        loci = self.path("loci.tsv")
        les = self.path("lesions.tsv")
        hits = self.path("hits.tsv")
        argv = ["parse", "--tblastn", self.hsp_file(), "--out-loci", loci, "--out-lesions", les, "--out-hits", hits] + list(extra)
        self.assertEqual(ts.main(argv), 0)
        return loci, les, hits

    def test_outputs_and_gene_inference(self):
        loci, les, hits = self.run_parse()
        h, rows = tsv_rows(read_text(loci))
        self.assertEqual(h, ts.LOCI_HEADER)
        self.assertEqual(len(rows), 4)
        # khong manifest: gene = qseqid, role = target
        self.assertTrue(all(r["gene"] == "ENAM_test" and r["role"] == "target" for r in rows))
        lesion = [r for r in rows if r["sseqid"] == "scaf_lesion"][0]
        self.assertEqual((lesion["n_stop_supported"], lesion["n_frameshift"]), ("1", "1"))
        lh, lrows = tsv_rows(read_text(les))
        self.assertEqual(lh, ts.LESIONS_HEADER)
        self.assertTrue(any(r["type"] == "frameshift" for r in lrows))
        hh, hrows = tsv_rows(read_text(hits))
        self.assertEqual(len(hrows), 7)
        self.assertNotIn("sseq", hh)

    def test_manifest_maps_gene_and_role(self):
        m = self.path("manifest.tsv")
        write_text(m, MANIFEST)
        loci, _, _ = self.run_parse(["--manifest", m])
        _, rows = tsv_rows(read_text(loci))
        self.assertTrue(all(r["gene"] == "ENAM" and r["accession"] == "XP_TEST.1" for r in rows))

    def test_split_qseqid_pipe_format(self):
        self.assertEqual(ts.split_qseqid("GULO|target|XP_1.1"), ("GULO", "target", "XP_1.1"))
        self.assertEqual(ts.split_qseqid("MEPE|control|XP_2.1"), ("MEPE", "control", "XP_2.1"))
        self.assertEqual(ts.split_qseqid("sp|P12345|ENAM_HUMAN"), ("sp|P12345|ENAM_HUMAN", "target", ""))


class TestExtract(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.prot, self.genome = make_synthetic_genome()
        self.fna = self.path("genome.fna.gz")
        with gzip.open(self.fna, "wt", encoding="utf-8", newline="") as f:
            for k, v in self.genome.items():
                ts.write_fasta_record(f, k, v)
        self.loci = self.path("loci.tsv")
        self.les = self.path("les.tsv")
        self.assertEqual(ts.main(["parse", "--tblastn", self.hsp_file(), "--out-loci", self.loci, "--out-lesions", self.les]), 0)

    def test_extract_pad_clip_and_strand(self):
        out = self.path("loci.fa")
        tr = self.path("tr.txt")
        rc = ts.main(["extract", "--fna", self.fna, "--loci", self.loci, "--pad", "10000", "--out", out, "--translate", tr])
        self.assertEqual(rc, 0)
        recs = {sid: seq for sid, _d, seq in ts.iter_fasta(io.StringIO(read_text(out)))}
        self.assertEqual(len(recs), 4)
        # pad lon hon scaffold -> cat ve [1, len]
        intact = [k for k in recs if k.startswith("ENAM_test|scaf_intact:")][0]
        self.assertTrue(intact.endswith(":1-{n}(+)".format(n=len(self.genome["scaf_intact"]))))
        self.assertEqual(recs[intact], self.genome["scaf_intact"].upper())
        # mach - duoc dao nguoc bo sung -> bang scaf_lesion (goc)
        rc_key = [k for k in recs if "scaf_lesion_rc" in k][0]
        self.assertTrue(rc_key.endswith("(-)"))
        self.assertEqual(recs[rc_key], self.genome["scaf_lesion"].upper())
        txt = read_text(tr)
        self.assertIn("## ENAM_test|scaf_lesion:", txt)
        self.assertIn("khung +1", txt)
        self.assertIn("khung +3", txt)

    def test_extract_small_pad_coordinates(self):
        out = self.path("loci.fa")
        self.assertEqual(ts.main(["extract", "--fna", self.fna, "--loci", self.loci, "--pad", "0", "--out", out]), 0)
        recs = {sid: seq for sid, _d, seq in ts.iter_fasta(io.StringIO(read_text(out)))}
        intact = [k for k in recs if "scaf_intact" in k][0]
        self.assertTrue(intact.endswith(":501-1040(+)"))
        self.assertEqual(recs[intact], self.genome["scaf_intact"][500:1040].upper())
        # vung trich doc dung khung +1 -> protein truy van (khong codon dung)
        self.assertEqual(ts.translate(recs[intact]), self.prot)

    def test_missing_scaffold_exit_1(self):
        bad = self.path("loci_bad.tsv")
        txt = read_text(self.loci).replace("scaf_intact", "scaf_missing")
        write_text(bad, txt)
        rc = ts.main(["extract", "--fna", self.fna, "--loci", bad, "--pad", "0", "--out", self.path("x.fa")])
        self.assertEqual(rc, 1)

    def test_translate_cli(self):
        fa = self.path("in.fa")
        write_text(fa, ">s1\nATGTAATTT\n")
        out = self.path("out.txt")
        self.assertEqual(ts.main(["translate", "--fasta", fa, "--out", out]), 0)
        self.assertIn("## s1", read_text(out))
        self.assertIn("1 codon dung tai nt 4", read_text(out))


LOCI_HDR = "\t".join(ts.LOCI_HEADER) + "\n"


def _loci_row(gene, role, sseqid, qcov, stops, fs, acc="XP_1.1", bits=100.0):
    vals = {
        "gene": gene, "role": role, "accession": acc, "qseqid": gene + "|" + role + "|" + acc,
        "sseqid": sseqid, "strand": "+", "s_start": 100, "s_end": 900, "n_hsp": 2, "qlen": 200,
        "qcov": qcov, "best_evalue": 1e-50, "max_pident": 60.0, "sum_bitscore": bits, "frames": "1,2",
        "n_stop_raw": stops, "n_stop_supported": stops, "n_frameshift": fs, "n_intron_like": 1, "n_in_frame_break": 0,
    }
    return "\t".join(ts._fmt_cell(vals[h]) for h in ts.LOCI_HEADER) + "\n"


class TestReport(TempDirMixin, unittest.TestCase):
    def run_report(self, loci_txt, lesions_txt=None, manifest_txt=None):
        loci = self.path("loci.tsv")
        write_text(loci, loci_txt)
        les = self.path("les.tsv")
        write_text(les, lesions_txt if lesions_txt is not None else "\t".join(ts.LESIONS_HEADER) + "\n")
        argv = ["report", "--loci", loci, "--lesions", les, "--out", self.path("r.md"), "--date", "2026-10-10"]
        if manifest_txt is not None:
            m = self.path("m.tsv")
            write_text(m, manifest_txt)
            argv += ["--manifest", m]
        self.assertEqual(ts.main(argv), 0)
        return read_text(self.path("r.md"))

    def test_control_pass_and_target_lesions(self):
        txt = self.run_report(LOCI_HDR + _loci_row("MEPE", "control", "s1", 0.95, 0, 0) + _loci_row("ENAM", "target", "s2", 0.6, 2, 1))
        self.assertIn("ĐẠT —", txt)
        self.assertNotIn("KHÔNG ĐẠT", txt)
        self.assertIn("2 codon dừng có chỗ dựa + 1 lệch khung", txt)
        self.assertIn("STRONG_INFERENCE", txt)
        self.assertIn("Thứ tự khả năng giải thích", txt)
        self.assertIn("Lỗi lắp ráp", txt)
        self.assertIn("Meredith 2014", txt)
        self.assertIn("\"Không tìm thấy\" ≠ \"không tồn tại\"", txt)
        self.assertIn("Cấm viết \"đã chứng minh\"", txt)

    def test_control_fail_blocks_interpretation(self):
        txt = self.run_report(LOCI_HDR + _loci_row("MEPE", "control", "s1", 0.95, 1, 0) + _loci_row("ENAM", "target", "s2", 0.9, 0, 0))
        self.assertIn("KHÔNG ĐẠT", txt)
        self.assertIn("không diễn giải được", txt)
        self.assertIn("chỉ để gỡ lỗi", txt)

    def test_no_control_blocks_interpretation(self):
        txt = self.run_report(LOCI_HDR + _loci_row("ENAM", "target", "s2", 0.9, 0, 0))
        self.assertIn("Không có gene đối chứng nào", txt)
        self.assertIn("không thấy tổn thương trong vùng căn được", txt)

    def test_manifest_missing_and_not_found(self):
        manifest = MANIFEST
        txt = self.run_report(LOCI_HDR + _loci_row("MEPE", "control", "s1", 0.95, 0, 0), manifest_txt=manifest)
        # ENAM co truy van nhung khong co locus -> "khong tim thay"; AMBN khong co truy van
        self.assertIn("| ENAM | XP_TEST.1 | 0 |", txt)
        self.assertIn("**không tìm thấy**", txt)
        self.assertIn("| AMBN |", txt)
        self.assertIn("truy vấn không có ở ngoài nhóm (not_found)", txt)

    def test_supported_lesions_table(self):
        les = "\t".join(ts.LESIONS_HEADER) + "\n" + "\t".join(
            ["ENAM", "target", "ENAM|target|XP_1.1", "s2", "+", "stop_codon", "91", "571", "1", "1", "0", "1", "1", "291", "query aa 91"]
        ) + "\n" + "\t".join(
            ["ENAM", "target", "ENAM|target|XP_1.1", "s2", "+", "stop_codon", "137", "700", "3", "3", "0", "0", "0", "88", "spurious"]
        ) + "\n"
        txt = self.run_report(LOCI_HDR + _loci_row("MEPE", "control", "s1", 0.95, 0, 0) + _loci_row("ENAM", "target", "s2", 0.9, 1, 0), lesions_txt=les)
        self.assertIn("| ENAM | s2 | + | stop_codon | 91 | 571 | 1→1 | query aa 91 | _điền tay_ |", txt)
        self.assertNotIn("spurious", txt)


# ---------------------------------------------------------------------------
# fetch-queries (mang gia lap)
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


def _fake_urlopen_factory(calls):
    def fake(req, timeout=60):
        url = req.full_url
        calls.append(url)
        q = dict(urllib_parse_qs(url))
        if "esearch.fcgi" in url:
            term = q["term"]
            if term.startswith("ENAM"):
                ids = ["111", "222"]
            elif term.startswith("MEPE"):
                ids = ["333"]
            else:
                ids = []
            return _FakeResp(json.dumps({"esearchresult": {"idlist": ids}}))
        if "efetch.fcgi" in url:
            ids = q["id"].split(",")
            fa = ""
            if "111" in ids:
                fa += ">XP_111.1 enamelin isoform X1 [Alligator]\nMKTAYIAKQR\n"
            if "222" in ids:
                fa += ">XP_222.1 enamelin isoform X2 [Alligator]\nMKTAYIAKQRQISFVKSHFSRQ\n"
            if "333" in ids:
                fa += ">XP_333.1 MEPE [Alligator]\nMSSSEE\n"
            return _FakeResp(fa)
        raise AssertionError("url la: " + url)

    return fake


def urllib_parse_qs(url):
    import urllib.parse as up

    qs = up.urlparse(url).query
    return {k: v[0] for k, v in up.parse_qs(qs).items()}


REGISTRY_HEADER = "dataset_id,name,source_url,version_or_release,downloaded_on,local_path,size_bytes,sha256,license,notes\n"
REGISTRY_ROW = "D01-x,old,http://x,v1,2026-09-07,data/x,1,abc,pub,note\n"


class TestFetchQueries(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.calls = []
        self._iv = ts.REQUEST_INTERVAL
        ts.REQUEST_INTERVAL = 0.0
        self.addCleanup(setattr, ts, "REQUEST_INTERVAL", self._iv)
        self._orig = ts._URLOPEN
        ts._URLOPEN = _fake_urlopen_factory(self.calls)
        self.addCleanup(setattr, ts, "_URLOPEN", self._orig)

    def test_fetch_picks_longest_and_records_not_found(self):
        out = self.path("q.faa")
        man = self.path("m.tsv")
        reg = self.path("registry.csv")
        write_text(reg, REGISTRY_HEADER + REGISTRY_ROW)
        rc = ts.main([
            "fetch-queries", "--genes", "ENAM,AMBN", "--controls", "MEPE", "--organism", "Alligator mississippiensis",
            "--email", "test@example.org", "--out", out, "--manifest", man, "--registry", reg, "--date", "2026-10-10",
        ])
        self.assertEqual(rc, 0)
        recs = list(ts.iter_fasta(io.StringIO(read_text(out))))
        ids = [r[0] for r in recs]
        self.assertEqual(ids, ["ENAM|target|XP_222.1", "MEPE|control|XP_333.1"])  # isoform dai nhat
        _, rows = tsv_rows(read_text(man))
        by_gene = {r["gene"]: r for r in rows}
        self.assertEqual(by_gene["ENAM"]["status"], "ok")
        self.assertEqual(by_gene["ENAM"]["n_candidates"], "2")
        self.assertEqual(by_gene["ENAM"]["length_aa"], "22")
        self.assertEqual(by_gene["AMBN"]["status"], "not_found")
        self.assertEqual(by_gene["MEPE"]["role"], "control")
        self.assertIn("refseq[Filter]", by_gene["ENAM"]["search_term"])
        # moi URL deu co email + tool
        self.assertTrue(all("email=test%40example.org" in u for u in self.calls))
        # registry: them dung 1 dong, dong cu giu nguyen, sha256 khop file
        reg_txt = read_text(reg)
        lines = reg_txt.splitlines()
        self.assertEqual(lines[0], REGISTRY_HEADER.strip())
        self.assertEqual(lines[1], REGISTRY_ROW.strip())
        self.assertEqual(len(lines), 3)
        sha = ts.sha256_file(out)
        self.assertIn(sha, lines[2])
        self.assertTrue(lines[2].startswith("D09-b5-queries-" + sha[:8]))
        self.assertIn("XP_222.1;XP_333.1", lines[2])
        # chay lai -> idempotent, khong them dong
        rc = ts.main([
            "fetch-queries", "--genes", "ENAM,AMBN", "--controls", "MEPE", "--organism", "Alligator mississippiensis",
            "--email", "test@example.org", "--out", out, "--manifest", man, "--registry", reg, "--date", "2026-10-10",
        ])
        self.assertEqual(rc, 0)
        self.assertEqual(len(read_text(reg).splitlines()), 3)

    def test_requires_email(self):
        with self.assertRaises(SystemExit):
            ts.main(["fetch-queries", "--email", "khong-phai-email", "--out", self.path("a"), "--manifest", self.path("b")])


# ---------------------------------------------------------------------------
# census: dem vung genome theo ho gene paralog
# ---------------------------------------------------------------------------

def _hsp(q, s, qstart, qend, sstart, send, bits, pident=60.0, qlen=300, frame=1):
    cols = [q, s, pident, qend - qstart + 1, 0, 0, qstart, qend, sstart, send, 1e-20, bits, qlen, 100000, frame, "A" * (qend - qstart + 1), "A" * (qend - qstart + 1)]
    return "\t".join(str(c) for c in cols)


class TestCensus(TempDirMixin, unittest.TestCase):
    def setUp(self):
        super().setUp()
        rows = [
            # vung 1 (scaf1 +, 1000-4500): HOXA1 thang ro (2 HSP, 2 exon cach 2 kb), HOXB1 chi homeobox
            _hsp("HOXA1|target|NP_1", "scaf1", 1, 120, 1000, 1360, 150.0, 70.0),
            _hsp("HOXA1|target|NP_1", "scaf1", 121, 300, 3400, 3940, 200.0, 80.0),
            _hsp("HOXB1|target|NP_2", "scaf1", 150, 210, 3500, 3683, 95.0, 65.0),
            # vung 2 (scaf1 +, 9000-9200): HOXA1 va HOXB1 gan bang nhau -> mo ho; cach vung 1 > 3 kb
            _hsp("HOXA1|target|NP_1", "scaf1", 150, 210, 9000, 9183, 90.0),
            _hsp("HOXB1|target|NP_2", "scaf1", 150, 210, 9001, 9184, 88.0),
            # vung 3 (scaf2 -): chi HOXB1, yeu (bits 50) -> khong "manh"
            _hsp("HOXB1|target|NP_2", "scaf2", 150, 210, 5000, 4817, 50.0, 55.0, 300, -1),
        ]
        self.tbl = self.path("t.tsv")
        write_text(self.tbl, "\n".join(rows) + "\n")

    def test_census_assigns_best_gene_and_flags_ambiguous(self):
        hsps = ts.read_hsps(self.tbl)
        rows = ts.census_loci(hsps, max_gap=3000)
        self.assertEqual(len(rows), 3)
        by = {r["locus_id"]: r for r in rows}
        v1 = by["scaf1:1000-3940(+)"]
        self.assertEqual(v1["best_gene"], "HOXA1")
        self.assertEqual(v1["best_bits"], 350.0)
        self.assertEqual(v1["second_gene"], "HOXB1")
        self.assertEqual(v1["n_queries_hit"], 2)
        self.assertEqual(v1["ambiguous"], 0)
        self.assertAlmostEqual(v1["best_qcov"], 1.0)
        v2 = by["scaf1:9000-9184(+)"]
        self.assertEqual(v2["best_gene"], "HOXA1")
        self.assertEqual(v2["ambiguous"], 1)  # 90 vs 88: bien < 10%
        v3 = by["scaf2:4817-5000(-)"]
        self.assertEqual(v3["best_gene"], "HOXB1")
        self.assertEqual(v3["second_gene"], "")
        self.assertEqual(v3["n_queries_hit"], 1)
        # sap theo best_bits giam dan
        self.assertEqual([r["best_gene"] for r in rows][0], "HOXA1")

    def test_census_cli_summary(self):
        out = self.path("census.tsv")
        summ = self.path("summary.tsv")
        rc = ts.main(["census", "--tblastn", self.tbl, "--out", out, "--summary", summ, "--min-bits", "80", "--min-qcov", "0.3"])
        self.assertEqual(rc, 0)
        h, rows = tsv_rows(read_text(out))
        self.assertEqual(h, ts.CENSUS_HEADER)
        self.assertEqual(len(rows), 3)
        _, srows = tsv_rows(read_text(summ))
        by = {r["gene"]: r for r in srows}
        # HOXA1: 2 vung thang (vung 1 manh + vung 2 mo ho, qcov 61/300 = 0.2 < 0.3 -> khong manh)
        self.assertEqual(by["HOXA1"]["n_loci"], "2")
        self.assertEqual(by["HOXA1"]["n_strong"], "1")
        self.assertEqual(by["HOXA1"]["n_unambiguous_strong"], "1")
        # HOXB1: chi vung 3 (bits 50 < 80) -> 1 locus, 0 manh
        self.assertEqual(by["HOXB1"]["n_loci"], "1")
        self.assertEqual(by["HOXB1"]["n_strong"], "0")

    def test_census_gap_splits_neighbouring_genes(self):
        hsps = ts.read_hsps(self.tbl)
        rows = ts.census_loci(hsps, max_gap=10000)  # gap lon -> vung 1 va 2 nhap lam mot
        self.assertEqual(len([r for r in rows if r["sseqid"] == "scaf1"]), 1)


# ---------------------------------------------------------------------------
# --help
# ---------------------------------------------------------------------------


class TestHelp(unittest.TestCase):
    def test_help_exit_zero(self):
        for argv in (["--help"], ["fetch-queries", "--help"], ["parse", "--help"], ["extract", "--help"], ["translate", "--help"], ["report", "--help"], ["census", "--help"]):
            with self.assertRaises(SystemExit) as cm:
                ts.main(argv)
            self.assertEqual(cm.exception.code, 0, argv)


# ---------------------------------------------------------------------------
# Tich hop voi BLAST+ that (bo qua neu khong co tblastn/makeblastdb)
# ---------------------------------------------------------------------------


@unittest.skipUnless(shutil.which("tblastn") and shutil.which("makeblastdb"), "can tblastn + makeblastdb tren PATH")
class TestIntegrationBlast(TempDirMixin, unittest.TestCase):
    def test_real_tblastn_roundtrip(self):
        prot, genome = make_synthetic_genome()
        fna = self.path("genome.fna")
        with open(fna, "w", newline="", encoding="utf-8") as f:
            for k, v in genome.items():
                ts.write_fasta_record(f, k, v)
        faa = self.path("q.faa")
        write_text(faa, ">ENAM|target|XP_TEST.1 test\n" + prot + "\n")
        db = self.path("db")
        subprocess.run(["makeblastdb", "-in", fna, "-dbtype", "nucl", "-out", db], check=True, stdout=subprocess.DEVNULL)
        tbl = self.path("t.tsv")
        with open(tbl, "w") as out:
            subprocess.run(["tblastn", "-query", faa, "-db", db, "-evalue", "1e-5", "-seg", "no", "-outfmt", ts.OUTFMT], check=True, stdout=out)
        loci = self.path("loci.tsv")
        les = self.path("les.tsv")
        self.assertEqual(ts.main(["parse", "--tblastn", tbl, "--out-loci", loci, "--out-lesions", les]), 0)
        _, rows = tsv_rows(read_text(loci))
        by = {r["sseqid"]: r for r in rows}
        self.assertEqual(set(by), {"scaf_intact", "scaf_lesion", "scaf_lesion_rc", "scaf_intron"})
        self.assertEqual((by["scaf_intact"]["n_stop_supported"], by["scaf_intact"]["n_frameshift"]), ("0", "0"))
        self.assertEqual((by["scaf_lesion"]["n_stop_supported"], by["scaf_lesion"]["n_frameshift"]), ("1", "1"))
        self.assertEqual((by["scaf_lesion_rc"]["n_stop_supported"], by["scaf_lesion_rc"]["n_frameshift"]), ("1", "1"))
        self.assertEqual(by["scaf_intron"]["n_stop_supported"], "0")
        self.assertGreaterEqual(int(by["scaf_intron"]["n_intron_like"]), 1)
        self.assertTrue(all(r["gene"] == "ENAM" and r["role"] == "target" for r in rows))


if __name__ == "__main__":
    unittest.main()
