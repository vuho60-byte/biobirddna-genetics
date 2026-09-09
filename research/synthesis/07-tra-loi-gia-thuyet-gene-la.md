# 07 — Trả lời giả thuyết "bộ gene lạ trong chim" bằng dữ liệu

Nguồn: `research/raw/R10-gene-dac-huu-chim.md` (26 claim, 21 C / 4 S / 3 U) và một lượt nghiên cứu sâu Perplexity (2026-09-09). Nhãn giữ nguyên của nguồn.

## 1. Có DNA "ngoại lai" trong bộ gene chim không? — CÓ, và đã biết nguồn

Chuyển gene ngang ở chim là hiện tượng **có thật**, nhưng con đường là **retrovirus**, không phải vi khuẩn và không phải nguồn ngoài Trái Đất:
- Gene envelope có nguồn gốc **thú** được tìm thấy trong retrovirus nội sinh của **ít nhất 15 trên 101 bộ gene chim** (Chen & Cui 2019 — S, một nguồn).
- Không tìm được ca chuyển gene ngang từ vi khuẩn sang chim.

Đây chính là loại "gene lạ nằm chung cho lớp chim" mà giả thuyết dự đoán, và nó có nguồn gốc đã truy được. Kết quả quét thành phần của chính dự án (`data/anomaly/`) cũng phù hợp: 0 trên 14.003 vùng mã hóa lệch cả ba thước đo.

## 2. Có allele trội gây yếu và gây chết không? — CÓ, ví dụ chuẩn ở gà

- **Creeper** (chân ngắn ở gà): allele **trội**, gây **chết ở thể đồng hợp**, tồn tại trong quần thể vì thể dị hợp sống được (C, xác nhận chéo ba lần, cơ chế phân tử đã giải).
- **Crest** ở gà: gây thoát vị não (C, hai bài cùng nhóm tác giả).
- **Araucana Rumpless** và **Ear-tufted** (C).

Đây đúng mô tả trong giả thuyết: gene khi trội thì gây đột biến và làm con vật yếu. Nhưng cơ chế đã biết: đột biến điều hòa hoặc mất chức năng ở gene phát triển, duy trì trong quần thể do **chọn lọc của người nuôi** giữ lại kiểu hình lạ, chứ không phải dấu vết của nguồn gốc khác.

## 3. Có nhóm gene chỉ riêng chim không? — có ở cấp thấp hơn, chưa ai đếm ở cấp lớp

- **Không tồn tại** công trình đếm gene đặc hữu toàn lớp Aves đối xứng với con số 274 gene **mất** (Lovell 2014). Đây là **khoảng trống nghiên cứu thật**, kiểm bằng hai đường độc lập.
- Có ở cấp thấp hơn: **10 gene chỉ có ở chim hót** (Wirthlin 2014) — liên quan trực tiếp tới chào mào.
- Gene mồ côi: chưa có tỉ lệ riêng cho chim; khung chung liên ngành là 10–30% (C, không riêng chim). Công cụ ORFanID có thể tự chạy cho *P. jocosus*.

## 4. Khoảng trống nghiên cứu cho họ Chào mào — xác nhận độc lập

Nghiên cứu sâu xác nhận: **chưa có** bản đồ vùng bảo tồn không mã hóa, **chưa có** track bảo tồn riêng, **chưa có** tái dựng trình tự tổ tiên, **chưa có** nghiên cứu gene sắc tố cho họ Pycnonotidae hoặc chi *Pycnonotus*. Bộ gene `GCA_013400435.1` mới chỉ được dùng làm tham chiếu trong các phân tích rộng.

Nghĩa là bản đồ lõi của dự án này là bản đầu tiên cho họ Chào mào.

## 5. Phương pháp chuẩn cho ba việc tiếp theo

| Việc | Công cụ chuẩn | Ràng buộc bắt buộc |
|---|---|---|
| Vùng tăng tốc riêng nhánh chào mào | **PhyloAcc**, bản **PhyloAcc-GT** xử lý bất tương thích cây gene và cây loài | mô hình trung tính từ vị trí thoái hóa bậc bốn; **phastBias** để loại chuyển đổi gene thiên vị GC; hiệu chỉnh đa kiểm định |
| Kiểm gene giàu lõi hơn mong đợi | hoán vị **giữ nguyên GC và bối cảnh bộ gene** | không dùng gán gene gần nhất; dùng cách định vùng điều hòa kiểu **GREAT** |
| Tìm mảnh vỡ gene men răng | **TBLASTN** với ortholog thú làm truy vấn | kiểm thủ công dịch khung và codon dừng sớm |

Riêng mục cuối: chưa ai tìm riêng ở bộ Sẻ, nên làm trên chào mào là việc chưa có tiền lệ.

## 6. Ghi chú độ tin cậy

Trong lượt tra cứu này, Perplexity trả **sai** số gene mất của Lovell 2014 (đưa ra 1.559 thay vì 274). Subagent tự đọc lại tóm tắt bài gốc và bỏ số sai. Đây là lần thứ hai một mô hình ngoài đưa số hoặc nguồn sai trong dự án; quy tắc kiểm chéo giữ nguyên.
