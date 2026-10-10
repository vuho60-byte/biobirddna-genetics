# p07_tooth_genes — việc B5 (brief X9): quét mảnh vỡ gene răng và GULO trên bộ gene chào mào

Bối cảnh đầy đủ: `research/briefs/X9-quet-gene-rang-chao-mao.md`; câu hỏi khoa học ở `docs/00` Q4 và
`research/synthesis/11-khung-co-the-chung-cung-to-tien-hay-cung-thiet-ke.md` mục 7. Phong cách CLI/TSV theo
`pipeline/p04_anomaly/composition_scan.py`; cơ chế resume `.params` sao chép nguyên văn từ
`pipeline/p02_core_map/run_a2.sh` (không import chéo — mỗi lối trong `pipeline/` tự chứa).

## Câu hỏi khoa học và logic phép thử

Meredith 2014 (R12 dòng #23, CONFIRMED): 6 gene men/ngà răng **ENAM, AMBN, AMELX, AMTN, DSPP, MMP20** mang
**đột biến bất hoạt dùng chung** ở hầu hết bộ chim đã khảo sát — khớp một lần mất răng khoáng hóa ở tổ tiên
chung chim hiện đại (~116 Ma). Chưa ai tìm riêng ở bộ Sẻ (tổng hợp 07 mục 5). Phép thử trên chào mào:

| Kết quả | Đọc |
|---|---|
| Có mảnh, tổn thương **khớp** vị trí Meredith 2014 | ủng hộ H1/H2 (cùng tổ tiên) trực tiếp trên loài của dự án |
| Có mảnh, tổn thương **không khớp** / riêng chào mào | vẫn tương thích H1 (giả gene tích thêm đột biến sau khi hỏng); **không** phải bằng chứng chống H1 |
| **Không tìm thấy** | **không kết luận**: lắp ráp mức scaffold (N50 ~218 kb), giả gene 116 Ma phân kỳ mạnh, ngưỡng e-value |
| Đối chứng dương **không đạt** | toàn bộ dòng gene mục tiêu **không diễn giải được** trong lần chạy đó |

Đối chứng dương mặc định: **MEPE** (cùng họ SCPP với gene răng, còn chức năng ở gà — docs/02 mục E) và
**SPARC** (gene khung xương bảo tồn cao). Chúng phải được tìm thấy với độ phủ ≥ 0,7, không codon dừng, không
lệch khung; nếu không, quy trình (truy vấn, tham số, bộ gene) có vấn đề trước khi nói gì về gene răng.

**GULO** đi kèm để lấp khoảng trống #19 của R15 (chim chỉ có dữ liệu enzyme 1969; trạng thái gene ở sẻ chưa
đọc từ trình tự). Dự đoán H1: nếu GULO hỏng ở chào mào, kiểu hỏng phải lồng được với sẻ họ hàng gần và khác
hẳn linh trưởng/chuột lang.

## Chạy

### Wrapper (khuyến nghị)

```
# Git Bash Windows + Docker Desktop (như run_a2.sh / run_last.sh):
bash pipeline/p07_tooth_genes/run_tooth_scan.sh --email ban@example.org
# WSL2 / Linux có tblastn trên PATH (conda env env/environment.yml đã thêm blast):
bash pipeline/p07_tooth_genes/run_tooth_scan.sh --mode local --email ban@example.org
# Đã có FASTA truy vấn (bỏ bước tải):
bash pipeline/p07_tooth_genes/run_tooth_scan.sh --queries my_queries.faa --manifest my_manifest.tsv
```

Thứ tự 8 bước, mỗi bước resume theo `<output>.params` chứa danh tính input (`--force` chạy lại hết):

```
1 fetch-queries  (python, host, CẦN MẠNG NCBI)  -> data/b5/queries.faa, queries_manifest.tsv, +1 dòng data/registry.csv
2 gzip -dc       bộ gene D01                    -> data/b5/genome.fna   (~1,2 GB — KHÔNG commit; makeblastdb không đọc .gz)
3 makeblastdb    (docker | local)               -> data/b5/blastdb/bulbul.*
4 tblastn        -seg no -max_target_seqs 50     -> data/b5/tblastn.tsv  (-outfmt cố định, xem dưới)
5 parse          (python)                        -> loci.tsv, lesions.tsv, hits.tsv
6 extract        (python, đọc .fna.gz theo luồng) -> loci.fa, loci_translate.txt   (để KIỂM THỦ CÔNG)
7 report         (python)                        -> report.md  (có phần cảnh báo bắt buộc)
8 sha256sum      các đầu ra nhỏ                  -> SHA256SUMS.txt
```

- `-seg no`: protein men răng giàu P/Q/S lặp; bộ lọc SEG mặc định sẽ che đúng vùng cần tìm.
- `-max_target_seqs 50`: giả gene thường vỡ ra nhiều scaffold ở lắp ráp mức scaffold.
- Docker image mặc định `quay.io/biocontainers/blast:2.16.0--h66d330f_5`: tag **suy từ** bioconda repodata
  (gói `blast-2.16.0-h66d330f_5`; biocontainers đặt tag `<version>--<build>`), **chưa `docker pull` kiểm**
  (mạng container bị chặn ở máy viết code). docs/02 ghi "cần image BLAST, chưa duyệt" → chủ dự án duyệt image
  rồi pull; tag khác dùng `--image`. Đường dẫn trong container luôn `/data/b5/...` (`/data` = `--data-dir`),
  nên đổi `--data-dir` không làm lệch đường dẫn (lesson 2026-09-09).
- `docker run --rm` là cờ dọn container sau khi chạy, **không** phải lệnh xóa file (giống `run_last.sh`).
- Ước lượng: 9 protein × bộ gene 1,2 Gb, tblastn 6 khung — cỡ chục phút tới một giờ tùy luồng; RAM vừa phải
  (BLAST quét theo volume ≤ 1 GB). Chưa đo trên dữ liệu thật.

### Lệnh con `tooth_scan.py` (Python 3.11 stdlib)

```
python -B tooth_scan.py fetch-queries --genes ENAM,AMBN,AMELX,AMTN,DSPP,MMP20,GULO --controls MEPE,SPARC \
    --organism "Alligator mississippiensis" --email X --out queries.faa --manifest queries_manifest.tsv \
    [--registry data/registry.csv] [--api-key K] [--all-proteins]
python -B tooth_scan.py parse --tblastn tblastn.tsv [--manifest queries_manifest.tsv] \
    --out-loci loci.tsv --out-lesions lesions.tsv [--out-hits hits.tsv] [--max-evalue 1e-5] \
    [--max-locus-gap 100000] [--flank 10] [--flank-min-ident 0.5] [--max-q-gap 20] [--max-fs-gap 40] [--min-intron 40]
python -B tooth_scan.py extract --fna genome.fna.gz --loci loci.tsv --pad 2000 --out loci.fa [--translate loci_translate.txt] [--role target]
python -B tooth_scan.py translate --fasta loci.fa [--out -]
python -B tooth_scan.py report --loci loci.tsv --lesions lesions.tsv [--manifest ...] --out report.md \
    [--min-control-cov 0.7] [--min-target-cov 0.3]
```

Bảng tblastn phải tạo bằng đúng chuỗi (lệnh `parse` kiểm số cột, sai thì dừng và in lại chuỗi này):

```
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sframe qseq sseq"
```

### 1. `fetch-queries`

esearch `GENE[Gene Name] AND "loài"[Organism] AND refseq[Filter]` (bỏ `refseq[Filter]` với `--all-proteins`) →
efetch FASTA → giữ **isoform dài nhất** (hòa: accession nhỏ hơn theo chuỗi — kết quả ổn định). Header ghi
`GENE|ROLE|ACCESSION` (`ROLE` = `target` | `control`) để `parse`/`report` tự nhận gene và vai trò.
`queries_manifest.tsv`: `qseqid gene role accession organism length_aa n_candidates status title search_term`;
gene không có ở ngoài nhóm → `status=not_found` (không crash — thiếu ortholog là thông tin). `--registry`:
thêm **đúng 1 dòng** theo header hiện có của `data/registry.csv`, `dataset_id = D09-b5-queries-<8 ký tự sha256>`,
idempotent (đã có cùng id thì không thêm), không xóa/ghi đè dòng cũ. **Không có accession nào viết sẵn trong
mã** — mọi accession đến từ NCBI tại thời điểm chạy và được ghi lại.

### 2. `parse` — HSP → locus → tổn thương

- Lọc `evalue <= --max-evalue`; nhóm HSP theo (truy vấn, scaffold, mạch); sắp theo tọa độ subject; nối thành
  **locus** khi khoảng cách subject ≤ `--max-locus-gap` (100 kb; MMP20 ở người trải ~48 kb).
- `loci.tsv`: `gene role accession qseqid sseqid strand s_start s_end n_hsp qlen qcov best_evalue max_pident
  sum_bitscore frames n_stop_raw n_stop_supported n_frameshift n_intron_like n_in_frame_break`.
  `qcov` = hợp các đoạn truy vấn được căn / `qlen`.
- **Codon dừng**: ký tự `*` trong `sseq`. `n_stop_raw` đếm tất cả; `n_stop_supported` chỉ đếm khi (a) vị trí
  truy vấn của `*` thuộc HSP **chủ sở hữu** (HSP có bitscore cao nhất phủ vị trí đó) và (b) mỗi bên có
  ≥ `--flank` cột căn không gap với tỉ lệ khớp ≥ `--flank-min-ident`. Lý do (a): trên dữ liệu tổng hợp, HSP
  thứ hai sau điểm lệch khung mở rộng ngược 8 aa vào vùng HSP thứ nhất ở **khung sai** và sinh `*` giả; HSP mở
  rộng vào intron cũng sinh `*` giả. Cột `s_nt` trong `lesions.tsv` = tọa độ nt 1-based (mạch +) của base thấp
  nhất trong codon.
- **Cặp HSP kề nhau trên truy vấn** (cùng locus) được phân loại theo **khoảng cách subject** `gap_nt`
  (mạch −: tính theo hướng gene): `> --min-intron` → `intron_like` (đổi khung giữa hai exon là bình thường,
  **không** phải tổn thương); cùng khung và ngắn → `in_frame_break`; khác khung, `|q_gap| ≤ --max-q-gap`,
  `gap_nt ≤ --max-fs-gap` (có thể âm = chồng lan) → `frameshift` ứng viên; còn lại → `unclassified` (để người xem).
  Không dùng chồng lan trên truy vấn để loại ranh exon: trên dữ liệu tổng hợp hai HSP hai bên intron 400 nt
  chồng lan tới 50 aa trên truy vấn vì cả hai mở rộng vào intron.
- `lesions.tsv`: `gene role qseqid sseqid strand type qpos s_nt frame_from frame_to gap_nt owned supported
  hsp_bitscore note` — một dòng cho mỗi `stop_codon` (kể cả không có chỗ dựa, `supported=0`), `frameshift`,
  `intron_like`, `in_frame_break`, `unclassified`.
- `hits.tsv` (tùy chọn): mọi HSP qua lọc, bỏ hai cột `qseq`/`sseq`.

### 3–4. `extract`, `translate`

`extract` đọc bộ gene `.fna(.gz)` theo luồng (một scaffold trong bộ nhớ), cắt `[s_start - pad, s_end + pad]`
(kẹp về `[1, len]`), **đảo ngược bổ sung** nếu mạch − để chuỗi ghi ra theo hướng đọc gene; header
`GENE|scaffold:start-end(mạch) pad=.. qseqid=..`. `--translate` ghi 3 khung dịch với số và tọa độ codon dừng
(nt 1-based trong đoạn trích). Thoát mã 1 nếu có locus không tìm thấy scaffold. `translate` làm việc đó cho
FASTA bất kỳ. Đây là bước **kiểm thủ công dịch khung và codon dừng** mà tổng hợp 07 mục 5 yêu cầu — công cụ
không thay được bước này.

### 5. `report`

`report.md` gồm: (1) bảng **đối chứng dương** ĐẠT/KHÔNG ĐẠT (độ phủ ≥ `--min-control-cov`, 0 stop có chỗ dựa,
0 lệch khung; không có đối chứng = KHÔNG ĐẠT); (2) bảng gene mục tiêu — locus tốt nhất theo `sum_bitscore`, cột
"Số locus" cho biết gene vỡ ra bao nhiêu mảnh; gene có trong manifest nhưng không locus → "**không tìm thấy**
(≠ không tồn tại)"; gene `not_found` ở ngoài nhóm → "chưa thử được"; (3) bảng tổn thương `supported=1` với cột
**"Khớp Meredith 2014?"** để trống cho người điền; (4) bước tiếp theo; (5) **phần cảnh báo bắt buộc** — không
tắt được.

## Thứ tự khả năng giải thích (BẮT BUỘC — in trong `report.md`)

Cho mỗi codon dừng / lệch khung, chỉ chuyển mục sau khi loại mục trước:

1. **Lỗi lắp ráp hoặc lỗi gọi base** (bộ gene mức scaffold; thay đổi 1 nt có thể là lỗi consensus) — kiểm bằng
   loài sẻ họ hàng gần có bộ gene tốt hơn (docs/05 B6).
2. **Giả tạo do căn trình tự** (HSP mở rộng vào intron/đầu đoạn với độ giống thấp) — cột `supported` lọc một
   phần; **không thay được kiểm thủ công** trên `loci_translate.txt`.
3. **Ranh exon thật** (`intron_like`) — không phải tổn thương.
4. **Tổn thương bất hoạt thật** — và chỉ khi **khớp** vị trí ở các bộ chim khác (bảng phụ lục Meredith 2014,
   chưa có trong repo) mới được gọi "cùng kiểu hỏng".

Nhãn tối đa cho mọi phát biểu rút từ công cụ này: **STRONG_INFERENCE** (một phương pháp, một bộ gene). Cấm
"đã chứng minh". Không nâng nhãn.

## Giới hạn

- **Chưa chạy trên dữ liệu thật** (bộ gene chào mào không có ở máy viết code; mạng tới NCBI/quay.io bị chặn).
  Mọi con số trong test là dữ liệu tổng hợp. Thời gian/RAM thật chưa đo.
- Chỉ phát hiện tổn thương **trong vùng tblastn căn được**; vùng gene đã mục nát quá mức (mất tín hiệu protein)
  hoặc nằm ở khoảng trống lắp ráp thì không thấy — kết quả âm **không kết luận**.
- Lệch khung **trùng ranh exon** không phân biệt được với ranh exon bình thường (cùng dấu hiệu: đổi khung +
  khoảng cách subject lớn). Indel bội số 3 không đổi khung → `in_frame_break`, không được coi là tổn thương.
- `*` do lỗi consensus lắp ráp không phân biệt được với `*` thật bằng một bộ gene; cần loài thứ hai.
- `_owner_by_qpos` chọn chủ sở hữu theo bitscore từng vị trí — nếu hai HSP cùng khung chồng lan (lặp đoạn),
  `*` có thể bị gán sai HSP; hiếm, để người kiểm qua `lesions.tsv` (`owned`).
- "Locus tốt nhất" trong báo cáo chọn theo `sum_bitscore`; tổn thương ở các locus khác vẫn có trong mục 3 và
  `lesions.tsv` — không bỏ sót, nhưng cột "Đọc tạm" chỉ nói về locus tốt nhất.
- `fetch-queries` lấy **isoform dài nhất** của RefSeq theo tên gene; tên gene khác ở loài ngoài nhóm (ví dụ
  ký hiệu LOC...) sẽ ra `not_found` — khi đó dùng `--all-proteins` hoặc tự cung cấp `--queries`.
- Docker image chưa pull kiểm; chế độ `--mode local` đã smoke test với BLAST+ 2.12.0 (Ubuntu), chuỗi `-outfmt`
  đọc được ở phiên bản đó. Phiên bản khác: test tích hợp sẽ báo nếu cột đổi.

## Test

`tests/test_tooth_scan.py` — unittest thuần, dữ liệu nhúng. Bảng HSP nhúng là **đầu ra thật của tblastn
2.12.0+** trên bộ gene tổng hợp (seed cố định, hàm `make_synthetic_genome` trong test): scaffold nguyên vẹn /
codon dừng tại aa 91 + xóa 1 nt tại codon 141 (mạch +) / cùng locus đảo ngược (mạch −) / hai exon cách intron
400 nt. Chạy:

```
python -B -m unittest discover -s pipeline/p07_tooth_genes/tests -v
```

32 test (đo 2026-10-10): dịch mã 3 khung + đảo mạch; đọc bảng, lọc e-value, sai số cột → dừng; nối locus (gap
100 kb → 4 locus; gap 100 nt → tách hai exon); codon dừng đúng aa 91, nt 571 (+) / 567 (−); `*` giả của HSP yếu
**không** tính (`owned=0`); lệch khung 1→3 / −1→−3; ranh exon → `intron_like`, không `frameshift`, `*` giả
trong intron không tính; `_flank_supported`, `_union_len`; CLI `parse` có/không manifest; `extract` kẹp pad ở
biên scaffold, đảo mạch −, khung +1 của vùng trích dịch ra đúng protein truy vấn, scaffold thiếu → exit 1;
`translate`; `report` khóa diễn giải khi đối chứng không đạt / không có đối chứng, in đủ cảnh báo, bảng tổn
thương chỉ dòng `supported=1`; `fetch-queries` với mạng **giả lập** (chọn isoform dài nhất, `not_found`,
registry thêm 1 dòng, idempotent, dòng cũ nguyên); `--help` mọi lệnh thoát 0; **test tích hợp** chạy
makeblastdb + tblastn thật (bỏ qua nếu không có trên PATH).

### Smoke test wrapper (đã làm 2026-10-10, `--mode local`, BLAST+ 2.12.0)

- Chạy **không** truyền `--data-dir`: báo thiếu bộ gene đúng đường dẫn mặc định `<repo>/data/ncbi/...` và exit 4
  (lesson 2026-09-09: kiểm đường dẫn mặc định, không chỉ đường dẫn ghi đè).
- Chạy với `--data-dir` trỏ vào bản sao `data/` có bộ gene tổng hợp nén `.fna.gz` + `--queries/--manifest`:
  đủ 8 bước, `report.md` báo đối chứng KHÔNG ĐẠT (manifest giả có MEPE `not_found`) và khóa diễn giải đúng;
  chạy lần hai → 7 STEP SKIP; `data/registry.csv` bản sao không đổi (chế độ `--queries` không ghi registry).
- `bash -n` sạch; grep `rm|rmdir|del|Remove-Item` chỉ khớp cờ `docker run --rm`.
