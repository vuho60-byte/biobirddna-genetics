# X7 — Mô phỏng "pha nước đường": sau bao nhiêu đời thì dấu vết lai biến mất

Chỉ tạo/sửa trong `pipeline/p05_admixture/**` và MỘT dòng thêm vào `tasks/primordial-dna-handoff.md`. **Không có lệnh xóa nào** (rm/rmdir/del/Remove-Item). Python 3.11 **stdlib** (được dùng `random`, `math`, `statistics`, `argparse`, `csv`, `json`). Không numpy/scipy. Không chạy docker. Không tải dữ liệu.

## Câu hỏi khoa học

Chủ dự án mô tả mô hình này: một bộ gene "gốc" trộn với một bộ gene "bản địa" rồi loãng dần qua các đời lai ngược, giống pha nước đường, giống khoai lang Nhật qua F1-F4. Câu hỏi của anh ấy là có tái dựng lại quá trình đó được không.

Câu hỏi trả lời được bằng mô phỏng là: **sau bao nhiêu đời thì một sự kiện lai trở nên KHÔNG THỂ phát hiện, và cái gì quyết định ngưỡng đó?**

Đây là mô hình chuẩn của introgression. Hai đại lượng suy giảm theo hai cách khác nhau, và đó chính là mấu chốt:

- **Tỉ lệ tổ tiên** giảm theo cấp số nhân dưới lai ngược liên tiếp: 1/2, 1/4, 1/8... Nhưng nếu quần thể lai giao phối trong nội bộ thay vì lai ngược thì tỉ lệ **không đổi**, chỉ dao động do trôi dạt.
- **Chiều dài đoạn tổ tiên** giảm do tái tổ hợp cắt nhỏ dần. Đoạn dài nghĩa là sự kiện lai gần đây; đoạn ngắn nghĩa là xa xưa. Kỳ vọng chiều dài đoạn sau `g` đời khoảng `1/(g·r)` Morgan.

Chính đại lượng thứ hai là "đồng hồ", và cũng chính nó đặt ra giới hạn: khi đoạn ngắn tới mức lẫn vào nhiễu thì không phương pháp nào tìm ra nữa.

## Việc — `pipeline/p05_admixture/dilution_sim.py`, ba lệnh con

1. `simulate --generations N --mode backcross|isolated --init-fraction 0.5 --chrom-length-cm 100 --n-chrom 30 --replicates 100 --seed 42 --out tracts.tsv --stats sim_stats.tsv`

   Mô phỏng Markov đơn giản trên mỗi nhiễm sắc thể: theo dõi tập các đoạn có nguồn gốc "gốc" (danh sách khoảng trên trục cM). Mỗi đời: sinh số điểm tái tổ hợp theo Poisson với kỳ vọng `chrom_length_cm/100` giao chéo mỗi nhiễm sắc thể, chọn ngẫu nhiên một trong hai bản sao ở mỗi khoảng giữa hai điểm chéo, rồi (chế độ `backcross`) ghép với một giao tử thuần bản địa, hoặc (chế độ `isolated`) ghép với một cá thể khác lấy từ chính quần thể mô phỏng.

   Xuất `tracts.tsv`: `replicate generation chrom start_cm end_cm len_cm`.
   Xuất `sim_stats.tsv`: mỗi đời một dòng — `generation mean_ancestry_fraction sd_ancestry_fraction n_tracts mean_tract_cm median_tract_cm max_tract_cm`.

2. `detect --stats sim_stats.tsv --min-tract-cm X --min-fraction Y --out detect.tsv`

   Với mỗi đời, đánh giá sự kiện lai còn phát hiện được hay không theo hai ngưỡng do người dùng đặt: tỉ lệ tổ tiên còn trên `--min-fraction` VÀ có ít nhất một đoạn dài hơn `--min-tract-cm`. Xuất cột `generation fraction_ok tract_ok detectable`.

   Trong README ghi rõ: hai ngưỡng này là **tham số người dùng đặt, không phải hằng số tự nhiên**; giá trị thực tế phải lấy từ tài liệu phương pháp (brief R11 đang thu thập), và cho tới khi có số đó thì mọi kết luận từ `detect` chỉ mang tính minh họa.

3. `report --stats sim_stats.tsv --detect detect.tsv --out report.md`

   Bảng theo đời, kèm phần diễn giải bằng lời trả lời đúng câu hỏi: sự kiện lai xảy ra cách đây bao nhiêu đời thì còn thấy được, và hai chế độ `backcross` với `isolated` khác nhau ra sao. Phải in cảnh báo bắt buộc: mô hình này bỏ qua chọn lọc, bỏ qua thay đổi kích thước quần thể, bỏ qua biến thiên tốc độ tái tổ hợp dọc nhiễm sắc thể, và giả định bản đồ di truyền đồng nhất — nên nó cho **bậc độ lớn**, không cho con số dùng để công bố.

## Test (`pipeline/p05_admixture/tests/test_dilution_sim.py`, unittest, dữ liệu nhúng, seed cố định)

- Chế độ `backcross`: tỉ lệ tổ tiên trung bình sau `g` đời xấp xỉ `0,5^(g+1)`... **hãy tự kiểm chứng công thức đúng bằng chính mô phỏng thay vì tin dòng này**; test phải khẳng định tỉ lệ giảm đơn điệu và khớp cấp số nhân trong sai số thống kê với ≥200 replicate.
- Chế độ `isolated`: tỉ lệ tổ tiên trung bình **không** giảm có hệ thống qua các đời (chỉ dao động). Đây là phép thử phân biệt hai chế độ, phải có.
- Chiều dài đoạn trung bình giảm đơn điệu theo số đời ở CẢ HAI chế độ.
- Không có đoạn nào vượt ra ngoài `[0, chrom_length_cm]`; tổng chiều dài đoạn luôn bằng `ancestry_fraction × tổng chiều dài nhiễm sắc thể` trong sai số làm tròn.
- Cùng `--seed` cho kết quả giống hệt nhau qua hai lần chạy.
- `--generations 0` trả về đúng tỉ lệ khởi đầu.

## Tiêu chí xong

- `python -B -m unittest discover -s pipeline/p05_admixture/tests -v` sạch; `--help` mỗi lệnh con thoát 0.
- Chạy thật một lần: `simulate --generations 20 --mode backcross` và `--mode isolated`, ghi kết quả vào `pipeline/p05_admixture/example/`, kích thước tổng dưới 5 MB.
- `pipeline/p05_admixture/README.md`: cách chạy, ý nghĩa từng cột, công thức dùng kèm nguồn nếu bạn xác minh được, và mục **Giới hạn** liệt kê đầy đủ những gì mô hình bỏ qua.
- Thêm đúng 1 dòng vào `tasks/primordial-dna-handoff.md`.

## Trả về ≤10 dòng
file tạo · dòng tổng kết unittest · số đời mà mô phỏng cho thấy tỉ lệ tụt dưới 1% ở chế độ backcross · xác nhận không lệnh xóa · giả định/giới hạn.
