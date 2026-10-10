# B5 — quét mảnh vỡ gene răng (và GULO) trên bộ gene chào mào

Ngày tạo: 2026-10-10 · công cụ: `pipeline/p07_tooth_genes/tooth_scan.py report` · đầu vào: `data/b5/loci.tsv`, `data/b5/lesions.tsv`

Chữ ký mọi kết quả: **STRONG_INFERENCE tối đa** — xem mục cuối trước khi trích dẫn bất kỳ dòng nào.

## 1. Đối chứng dương

| Gene | Accession | Locus tốt nhất | Độ phủ | Stop (chỗ dựa) | Lệch khung | Kết quả |
|---|---|---|---|---|---|---|
| MEPE | XP_059578272.1 | — | 0 | — | — | KHÔNG ĐẠT (không tìm thấy) |
| SPARC | XP_006270448.1 | VWYP01029506.1:237260-239474(+) | 0.95 | 0 | 0 | ĐẠT |

**Kết luận đối chứng:** KHÔNG ĐẠT — dòng gene mục tiêu bên dưới **không diễn giải được** cho tới khi sửa.

## 2. Gene mục tiêu

| Gene | Accession truy vấn | Số locus | Locus tốt nhất | Độ phủ | %id max | e-value | Stop thô / có chỗ dựa | Lệch khung | Ranh exon | Đọc tạm |
|---|---|---|---|---|---|---|---|---|---|---|
| AMBN | XP_019346383.1 | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |
| AMELX | XP_014463267.1 | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |
| AMTN | XP_014451634.1 | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |
| DSPP | XP_059578265.1 | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |
| ENAM | XP_006257995.1 | 0 | — | 0 | — | — | — | — | — | **không tìm thấy** (≠ không tồn tại; xem mục cuối) |
| GULO | — | — | — | — | — | — | — | — | — | truy vấn không có ở ngoài nhóm (not_found) → chưa thử được |
| MMP20 | XP_006272555.1 | 11 | VWYP01031341.1:429909-554428(-) | 0.58 | 70.27 | 3.47e-21 | 13 / 2 | 0 | 11 | có mảnh, độ phủ đạt; 2 codon dừng có chỗ dựa + 0 lệch khung ứng viên → kiểm thủ công, rồi đối chiếu Meredith 2014 |

> Đối chứng dương không đạt: bảng trên **chỉ để gỡ lỗi**, không dùng cho bất kỳ phát biểu nào.

## 3. Tổn thương ứng viên (chỉ dòng `supported=1`; đầy đủ trong `data/b5/lesions.tsv`)

| Gene | Scaffold | Mạch | Loại | Vị trí aa truy vấn | Toạ độ nt | Khung | Ghi chú | Khớp Meredith 2014? |
|---|---|---|---|---|---|---|---|---|
| MMP20 | VWYP01031341.1 | - | stop_codon | 173 | 552950 | -2→-2 | query aa 173 | _điền tay_ |
| MMP20 | VWYP01031341.1 | - | stop_codon | 271 | 479973 | -1→-1 | query aa 271 | _điền tay_ |

Cột cuối để trống có chủ đích: so với bảng phụ lục Meredith 2014 là bước **người làm**, không tự động.

## 4. Bước tiếp theo (theo thứ tự)

1. `tooth_scan.py extract --fna <genome.fna.gz> --loci data/b5/loci.tsv --pad 2000 --translate loci_translate.txt` rồi đọc 3 khung dịch quanh mỗi tổn thương.
2. Lấy bảng phụ lục Meredith 2014 (DOI 10.1126/science.1254390) → đối chiếu vị trí từng tổn thương; chỉ dòng **khớp** mới ghi "cùng kiểu hỏng".
3. Chạy lại cùng truy vấn trên một bộ gene sẻ chất lượng cao (sẻ vằn / sẻ ngô lớn — docs/05 B6) để tách lỗi lắp ráp khỏi tổn thương chung.
4. Ghi kết quả vào `research/synthesis/` với nhãn; **không** sửa `research/raw/`.

## Cách đọc kết quả — BẮT BUỘC đọc trước khi dùng bất kỳ dòng nào ở trên

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
