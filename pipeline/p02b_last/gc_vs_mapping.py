"""Ti le can duoc theo ham luong GC cua phan tu bao ton ngan.

Vi sao co script nay: sau B1, 97,74% phan tu nhom 20-49 bp van khong can duoc
sang bo gene chao mao. Co ba cach giai thich (san phat hien / lo hong ban lap
rap / vang mat sinh hoc that) va truoc day chua phan dinh duoc.

Bornelov et al. 2017 Genome Biology (DOI 10.1186/s13059-017-1231-1, da kiem
qua PubMed) cho thay 137 gene ga tung bi tuyen bo "mat trong tien hoa" thuc
ra chi vang mat khoi BAN LAP RAP vi nam trong vung giau GC do phu kem. Neu
co che do cung tac dong len du an nay thi phan tu giau GC phai kho can hon
mot cach co he thong. Day la phep thu cho du doan do.

Dau ra: bang ti le can duoc theo bin GC, tach rieng theo nhom chieu dai.

Chay:
  python -B pipeline/p02b_last/gc_vs_mapping.py \
    --short-fa data/a2_last/short.fa \
    --mapped-bed data/a2_last/bulbul_short.bed \
    --out data/a2_last/gc_vs_mapping.tsv
"""

import argparse
import gzip
import sys

GC_BINS = [
    (0.00, 0.30), (0.30, 0.35), (0.35, 0.40), (0.40, 0.45), (0.45, 0.50),
    (0.50, 0.55), (0.55, 0.60), (0.60, 0.70), (0.70, 1.01),
]
LEN_BINS = [("20_49", 20, 49), ("50_99", 50, 99)]
MIN_BIN_COUNT = 500  # duoi nguong nay khong in ti le (nhieu qua lon)


def smart_open(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") \
        else open(path, encoding="utf-8")


def read_mapped_ids(bed_path):
    ids = set()
    with smart_open(bed_path) as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) >= 4:
                ids.add(parts[3])
    return ids


def iter_fasta(handle):
    header, chunks = None, []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def gc_fraction(seq):
    """GC tren so base XAC DINH (bo qua N va ky tu la). None neu khong co base nao."""
    upper = seq.upper()
    acgt = sum(1 for c in upper if c in "ACGT")
    if acgt == 0:
        return None
    return sum(1 for c in upper if c in "GC") / acgt


def gc_bin_index(value):
    for i, (lo, hi) in enumerate(GC_BINS):
        if lo <= value < hi:
            return i
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--short-fa", required=True, help="FASTA phan tu ngan da loc (short.fa)")
    ap.add_argument("--mapped-bed", required=True, help="BED phan tu da can duoc (cot 4 = ma phan tu)")
    ap.add_argument("--out", help="TSV ket qua; khong co thi in ra stdout")
    args = ap.parse_args(argv)

    mapped = read_mapped_ids(args.mapped_bed)

    # counts[(len_label, gc_bin_index)] = [so can duoc, tong]
    counts = {}
    skipped_len = 0
    skipped_gc = 0
    with smart_open(args.short_fa) as f:
        for name, seq in iter_fasta(f):
            length = len(seq)
            len_label = None
            for label, lo, hi in LEN_BINS:
                if lo <= length <= hi:
                    len_label = label
                    break
            if len_label is None:
                skipped_len += 1
                continue
            gc = gc_fraction(seq)
            if gc is None:
                skipped_gc += 1
                continue
            idx = gc_bin_index(gc)
            if idx is None:
                skipped_gc += 1
                continue
            key = (len_label, idx)
            cell = counts.setdefault(key, [0, 0])
            cell[1] += 1
            if name in mapped:
                cell[0] += 1

    out = open(args.out, "w", encoding="utf-8", newline="\n") if args.out else sys.stdout
    try:
        out.write("len_bin\tgc_lo\tgc_hi\ttotal\tmapped\tmapped_pct\n")
        for label, _, _ in LEN_BINS:
            for i, (lo, hi) in enumerate(GC_BINS):
                m, t = counts.get((label, i), [0, 0])
                if t < MIN_BIN_COUNT:
                    continue
                out.write("%s\t%.2f\t%.2f\t%d\t%d\t%.2f\n"
                          % (label, lo, hi, t, m, 100.0 * m / t))
        out.write("\nmetric\tvalue\n")
        out.write("skipped_out_of_len_range\t%d\n" % skipped_len)
        out.write("skipped_no_acgt\t%d\n" % skipped_gc)
        out.write("min_bin_count_printed\t%d\n" % MIN_BIN_COUNT)
    finally:
        if args.out:
            out.close()
            print("da ghi " + args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
