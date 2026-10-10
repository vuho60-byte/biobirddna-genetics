# p08_silent_sites — việc B11 (brief X10): phép thử vị trí câm, Hox người so với chào mào

Bối cảnh: `research/briefs/X10-phep-thu-vi-tri-cam-hox.md`; câu hỏi gốc ở `docs/00` Q4 và tổng hợp 11. Phong cách
CLI/TSV theo `pipeline/p07_tooth_genes/tooth_scan.py` (không import chéo — mỗi lối tự chứa).

## Phép thử, viết kỳ vọng trước khi nhìn số

Người và chim dùng chung cụm Hox và bộ gene công cụ phát triển. Hai cách giải thích cho hai dự đoán khác nhau tại
**vị trí câm**: base thứ 3 của codon thuộc họ 4-fold (CTN Leu, GTN Val, TCN Ser, CCN Pro, ACN Thr, GCN Ala, CGN Arg,
GGN Gly), nơi mọi thay đổi đều không đổi axit amin nên chọn lọc lên protein không giữ được nó.

| Kịch bản | Protein (vùng lõi) | `ident_4fold_core` | `jc_4fold_core` |
|---|---|---|---|
| Thừa kế từ tổ tiên chung ~320 Ma dưới chọn lọc (H1/H2) | giống cao | tiến về `chance_4fold_core` | lớn hoặc `inf` |
| Chép/ghép gần đây giữa hai bộ gene | bất kỳ | gần 1,0 | gần 0 |

Đối chứng: gene housekeeping phải cho cùng bức tranh với Hox. Chỉ khi Hox lệch khỏi đối chứng theo hướng "giống hệt
vị trí câm" mới có tín hiệu đáng xem.

## Chạy

```
# 1 (CẦN MẠNG NCBI): CDS người RefSeq cho 39 HOX + 11 gene công cụ + 10 gene đối chứng
python -B pipeline/p08_silent_sites/silent_sites.py fetch-human --email ban@example.org \
    --out data/b11/human_cds.fna --manifest data/b11/human_manifest.tsv --registry data/registry.csv
# 2: so với CDS chào mào (D01)
python -B pipeline/p08_silent_sites/silent_sites.py compare --human data/b11/human_cds.fna \
    --bulbul data/ncbi/GCA_013400435.1/GCA_013400435.1_ASM1340043v1_cds_from_genomic.fna.gz \
    --out data/b11/silent_sites.tsv --alignments data/b11/alignments.txt
# 3: báo cáo
python -B pipeline/p08_silent_sites/silent_sites.py report --tsv data/b11/silent_sites.tsv --out data/b11/report.md
```

Chỉ Python 3.11 stdlib. Không docker. Thời gian: căn toàn cục thuần Python cho ~60 cặp protein 300–500 aa, cỡ vài
chục giây tới vài phút. CDS chào mào đọc theo luồng, chỉ giữ gene cần.

### `fetch-human`
esearch nuccore `GENE[Gene Name] AND "Homo sapiens"[Organism] AND refseq[Filter] AND biomol_mrna[PROP]` → efetch
`rettype=fasta_cds_na` → chỉ nhận bản ghi có `[gene=GENE]` đúng tên, CDS hợp lệ (chia hết 3, ≥ 30 codon, không codon
dừng nội) → giữ **CDS dài nhất**. Header ghi `GENE|human|NM_acc|NP_acc [gene=GENE] [protein_id=..]`. Manifest:
`qseqid gene class organism mrna_acc protein_id length_nt n_candidates status title search_term`. `--registry` thêm 1
dòng `D10-b11-human-cds-<sha8>` theo header hiện có, idempotent, không xóa/ghi đè. Mặc định gọi 52 tên HOX (A1–D13);
người chỉ có 39 nên 13 tên sẽ `not_found` — đó là kỳ vọng, không phải lỗi. **Không có accession nào viết sẵn trong mã.**

### `compare`
- Khớp gene theo tên không phân biệt hoa/thường (`HOXA1` người ↔ `[gene=Hoxa1]` chào mào). Mỗi bên giữ CDS dài nhất
  (nhiều isoform / mô hình gene). Bỏ codon dừng cuối. `status` ghi lý do loại: `human_/bulbul_` + `len_not_multiple_of_3`
  / `too_short` / `internal_stop`; `no_bulbul_cds` khi chú giải chào mào không có tên gene đó.
- Dịch protein → căn toàn cục Needleman–Wunsch/Gotoh (BLOSUM62 nhúng trong mã, gap mở 11, kéo dài 1) → ánh xạ về
  codon (chỉ cột không gap ở cả hai bên).
- **Cột lõi**: cột nằm trong cửa sổ ±`--window` (10) cột có tỉ lệ aa giống ≥ `--min-ident` (0,5). Loại vùng căn mơ hồ
  (đuôi Hox ngoài homeodomain, vùng lặp poly-Ala/Pro), nơi codon có thể bị ghép sai cặp.
- `silent_sites.tsv`, mỗi gene một dòng: `gene class human_id bulbul_id human_len_nt bulbul_len_nt aln_len aln_score
  n_codon_all aa_ident_all ident_p3_all n_4fold_all ident_4fold_all n_codon_core aa_ident_core ident_p1_core
  ident_p2_core ident_p3_core n_4fold_core ident_4fold_core chance_4fold_core jc_4fold_core gc3_human gc3_bulbul status`.
  `chance_4fold_core` = Σ_b f_người(b)·f_chào mào(b) trên các vị trí 4-fold lõi (mức giống ngẫu nhiên theo thành
  phần thật, không phải 0,25). `jc_4fold_core` = −¾·ln(1 − 4/3·p), `inf` khi p ≥ 0,75 (bão hòa). `status=ok` cần
  ≥ 20 vị trí 4-fold lõi; ít hơn → `few_4fold_sites`, không xét.
- `--alignments`: căn protein hai dòng + dòng `core` (`#` lõi, `.` ngoài lõi, khoảng trắng = gap) để kiểm thủ công.
- Lớp gene: `^HOX[ABCD]\d+$` → `hox`; danh sách công cụ mặc định (PAX6, SOX9, RUNX2, SHH, BMP4, TBX5, TBX4, PITX1,
  FGF8, GLI3, HAND2) → `toolkit`; còn lại → `control`; `--classes gene<TAB>class` ghi đè.

### `report`
Bảng theo lớp với trung vị, mục "So Hox với đối chứng", verdict từng gene (`--hi` 0,9: "gần giống hệt → dấu hiệu
chép/ghép gần đây, kiểm lại"; trong ±`--band` 0,15 quanh ngẫu nhiên: "khớp thừa kế lâu đời dưới chọn lọc"; giữa hai
mức: "thừa kế + ràng buộc codon"; dưới ngẫu nhiên: "kiểm căn trình tự"), và **phần cảnh báo bắt buộc** không tắt được.

## Giới hạn (in trong `report.md`)
1. Một cặp loài, một bộ gene mức scaffold; CDS chào mào là mô hình gene dự đoán — exon sai hạ `aa_ident` và kéo
   `ident_4fold` về ngẫu nhiên giả tạo; vì vậy dùng cột lõi và ngưỡng ≥ 20 vị trí. Lặp với gà/sẻ vằn trước khi vượt
   STRONG_INFERENCE.
2. GC khác nhau giữa hai loài làm mức ngẫu nhiên lệch khỏi 0,25; so với `chance_4fold_core`.
3. Codon usage bias / yếu tố điều hòa trong CDS giữ vị trí câm trên ngẫu nhiên một chút; không tạo ra ~1,0 trên cả gene.
4. Không phân biệt "thừa kế" với một kịch bản thiết kế cố tình ngẫu nhiên hóa vị trí câm; kịch bản đó không còn là
   "ghép/chép" và phải tự nêu dự đoán riêng (ô H3 trong `docs/00`).
5. Căn toàn cục không mô hình intron/exon; vùng gap lớn (isoform khác) bị loại khỏi thống kê, không gây sai nhưng làm
   giảm `n_codon`.
6. **Chưa chạy dữ liệu thật** (NCBI bị chặn ở máy viết code; D01 không có ở đây). Mọi số trong test là tổng hợp.

## Test
`tests/test_silent_sites.py` — unittest thuần, dữ liệu nhúng, mạng giả lập. Chạy:

```
python -B -m unittest discover -s pipeline/p08_silent_sites/tests -v
```

20 test (2026-10-10): BLOSUM62 đối xứng, đường chéo dương, W/W = 11; căn giống hệt không gap, chèn 3 aa → 3 gap đúng
chiều; bảng mã; 4-fold (CTN đúng, TTR sai, Ile sai, ký tự mơ hồ sai); JC69 (0 → 0, ≥ 0,75 → inf, 0,3 → −¾ ln 0,6);
`strip_stop`/`check_cds`; CDS giống hệt → 1,0; đổi base 3 mọi codon 4-fold → protein giữ, `ident_4fold` 0, JC inf;
ngẫu nhiên hóa → `ident_4fold` < 0,6, chance trong (0,15; 0,6), JC > 0,3; mặt nạ lõi loại đúng khối 30 aa thay bằng
axit amin tương tự (cột 45–74 ngoài lõi, hai bên lõi); `codon_columns` bỏ cột gap; `gene_class`; `compare`
end-to-end (isoform dài nhất, `bulbul_internal_stop`, `no_bulbul_cds`, thứ tự lớp, file căn); `report` đủ verdict 4
loại và cảnh báo; `fetch-human` mạng giả lập (CDS dài nhất đúng gene, bỏ gene sai, `not_found`, registry 1 dòng
idempotent, dòng cũ nguyên); `--help` 4 lệnh thoát 0.
