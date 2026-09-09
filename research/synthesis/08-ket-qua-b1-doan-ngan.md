# 08 — Kết quả B1: căn các đoạn bảo tồn ngắn bằng LAST

Ngày chạy: 2026-09-09. Đầu vào: `data/a2/elements.{fa,bed}` (2.947.588 phần tử bảo tồn trên gà galGal4), bộ gene chào mào `GCA_013400435.1`. Kịch bản: `pipeline/p02b_last/run_last.sh --chunk-bases 45000000 --threads 6`. Nhật ký: `data/a2_last/run.log`.

## 1. Vì sao phải làm bước này

Vòng A2 dùng minimap2 và chỉ căn được 8,43% số phần tử. Tỉ lệ căn được tăng đều theo chiều dài phần tử: 0,25% ở nhóm 20-49 bp, 4,73% ở nhóm 50-99 bp, lên tới 74,03% ở nhóm từ 500 bp. Nhóm dưới 100 bp chiếm 2.311.288 phần tử, tức 78,4% tổng số. Nếu tỉ lệ thấp ở nhóm ngắn là do công cụ chứ không do sinh học thì mọi con số của dự án đang lệch một cách có hệ thống. B1 kiểm đúng giả thuyết đó bằng LAST, công cụ chuẩn cho trình tự ngắn giữa hai loài xa nhau. [CONFIRMED — số liệu từ `data/a2/map_stats.tsv` của chính dự án]

## 2. Tỉ lệ căn được: LAST so với minimap2

| Nhóm chiều dài | Số phần tử | minimap2 | LAST | Bội số |
|---|---|---|---|---|
| 20-49 bp | 1.770.154 | 0,25% | 2,26% | 9,0 |
| 50-99 bp | 541.134 | 4,73% | 15,72% | 3,3 |
| Cả hai nhóm | 2.311.288 | 1,30% | 5,41% | 4,2 |

[CONFIRMED — `data/a2/map_stats.tsv` và `data/a2_last/map_stats_short.tsv`]

LAST tìm thêm được 95.169 phần tử ngắn mà minimap2 bỏ sót. Phân bố mức giống nhau của các đoạn mới: 26 đoạn dưới 0,80, còn lại đều từ 0,80 trở lên, trong đó 24.366 đoạn giống hệt 100%.

## 3. Kiểm chéo độc lập giữa hai công cụ

29.827 phần tử ngắn được cả hai công cụ căn được. Đây là phép kiểm chéo thật vì hai công cụ dùng thuật toán tìm hạt giống khác hẳn nhau.

| Tiêu chí đối chiếu | Số phần tử | Tỉ lệ |
|---|---|---|
| Cùng scaffold | 29.729 | 99,67% |
| Tọa độ trùng khít từng base | 25.611 | 85,87% |
| Lệch không quá 50 bp và có chồng lấn | 29.696 | 99,56% |
| Khác scaffold hoàn toàn | 98 | 0,33% |

[CONFIRMED — tính trực tiếp từ `data/a2/bulbul_core.bed` và `data/a2_last/bulbul_short.bed`]

Đây là bằng chứng mạnh rằng cả hai công cụ đang tìm ra tương đồng thật chứ không phải khớp ngẫu nhiên. Đồng thời nó cũng cho thấy 0,33% số ca hai công cụ đặt phần tử lên hai scaffold khác nhau, phần này chưa giải quyết và nên loại khi dùng cho phân tích cần độ chính xác vị trí.

## 4. Ảnh hưởng lên các con số đã công bố

| Chỉ số | Chỉ đoạn dài | Đầy đủ | Thay đổi |
|---|---|---|---|
| Số phần tử trong lõi | 248.461 | 343.781 | +38,4% |
| Lõi sau gộp trùng (bp) | 61.498.037 | 66.627.572 | +8,3% |
| Tỉ lệ trên bộ gene | 6,00% | 6,50% | +0,50 điểm |
| Vùng mã hóa nằm trong lõi | 59,8% | 61,32% | +1,5 điểm |
| Mức giàu vùng mã hóa | 9,96 lần | 9,43 lần | −0,53 |

[CONFIRMED — `data/a2/overlap_stats.tsv` và `data/a2_last/overlap_stats_all.tsv`, cùng một chương trình `overlap_stats.py`, cùng cách tính giao từng base]

**Kết luận quan trọng nhất của B1: các con số cũ không bị lật.** Tôi đã lo rằng việc bỏ sót 78% số phần tử làm lệch mọi kết quả. Thực tế bản đồ lõi chỉ lớn thêm 8,3% về số base, và ba phát biểu chính vẫn đứng: lõi chiếm khoảng 6-6,5% bộ gene, hơn 60% vùng mã hóa nằm trong lõi, và vùng mã hóa giàu lõi hơn kỳ vọng khoảng 9-10 lần. Lý do là các phần tử ngắn tuy đông nhưng mỗi cái chỉ vài chục base, nên đóng góp ít vào tổng số base.

Mức giàu giảm nhẹ từ 9,96 xuống 9,43 vì các đoạn ngắn mới phân bố ra ngoài vùng mã hóa nhiều hơn phần cũ. Đây là chiều đúng với kỳ vọng: phần tử bảo tồn ngắn phần lớn là yếu tố điều hòa nằm ngoài gene.

## 5. Điều B1 KHÔNG giải quyết

97,74% số phần tử nhóm 20-49 bp vẫn không căn được, kể cả với công cụ chuyên dụng. Tỉ lệ căn được vẫn tăng đều theo chiều dài. [CONFIRMED — bảng mục 2]

Ba cách giải thích, chưa phân định được bằng dữ liệu hiện có:

1. **Sàn phát hiện.** Một đoạn 25 bp đã phân kỳ khoảng 90 triệu năm không còn đủ thông tin để bất kỳ thuật toán nào định vị chắc chắn. Đây là giới hạn lý thuyết, không phải lỗi công cụ. [STRONG_INFERENCE — chiều tăng đều theo chiều dài ở cả hai công cụ độc lập ủng hộ cách này]
2. **Lỗ hổng bản lắp ráp.** Bộ gene chào mào ở mức scaffold, độ liên tục thấp; phần tử rơi vào khoảng trống thì không thể căn. [UNVERIFIED — chưa đo tỉ lệ N và chưa đối chiếu với bản lắp ráp chất lượng cao của loài chim khác]
3. **Vắng mặt sinh học thật.** Các đoạn đó thực sự mất hoặc biến đổi ở dòng chào mào. [UNVERIFIED]

Phép thử phân định: chạy lại đúng quy trình này với một loài chim sẻ có bộ gene mức nhiễm sắc thể. Nếu tỉ lệ căn được ở nhóm 20-49 bp tăng vọt thì nguyên nhân là chất lượng bản lắp ráp; nếu vẫn thấp như nhau thì nguyên nhân là sàn phát hiện. Chưa làm.

B1 cũng không nói các đoạn này làm nhiệm vụ gì, không phân định ba giả thuyết nguồn gốc H1/H2/H3, và không sửa lỗi bảng xếp hạng gene bị thổi phồng do mô hình gene cụt.

## 6. Đầu ra

| Tệp | Nội dung |
|---|---|
| `data/a2_last/bulbul_short.bed` | 125.147 phần tử ngắn, tọa độ chào mào |
| `data/a2_last/bulbul_core_all.bed` | 343.781 phần tử, bản đồ lõi đầy đủ |
| `data/a2_last/map_stats_short.tsv` | Tỉ lệ căn theo nhóm chiều dài và theo mức giống nhau |
| `data/a2_last/overlap_stats_all.tsv` | Giao từng base với chú giải, tính trên bản đồ đầy đủ |
| `data/a2_last/aln.maf` | Kết quả căn thô của LAST, 54 MB |

Tài nguyên đã dùng: 23 mảnh bộ gene 45 Mbp, khoảng 35 giây mỗi mảnh, tổng 13 phút. Đỉnh bộ nhớ container quan sát được 206 MiB trên 1,9 GB khả dụng, không hề chạm trần như lo ngại ban đầu.

## 7. Một lỗi đã sửa trong lúc chạy

`run_last.sh` tính thư mục gốc dự án lùi một cấp trong khi nó nằm sâu hai cấp, nên `DATA_DIR` trỏ vào `pipeline/data` không tồn tại và kịch bản thoát mã 4 ngay lệnh kiểm tệp đầu vào. 14 test đơn vị không bắt được vì chúng chỉ kiểm phần Python. Lần chạy thử trên bộ dữ liệu nhỏ cũng không bắt được vì tôi truyền `--data-dir` tường minh. Đã sửa ở commit 64bd352 kèm bài học trong `tasks/lessons.md`.
