# 09 — Mô hình "pha nước đường": hai kịch bản có số phận trái ngược

Ngày: 2026-09-09. Nguồn: `research/raw/R11-introgression-ghost-lineage.md` (phương pháp) và `pipeline/p05_admixture/` (mô phỏng, 37/37 test đạt). Câu hỏi gốc do chủ dự án đặt: một bộ gene "gốc" trộn với bộ gene "bản địa" rồi loãng dần qua các đời — có tái dựng lại được quá trình đó không?

## 1. Mô hình có tên trong khoa học

Cái chủ dự án mô tả bằng hình ảnh pha nước đường chính là **introgression** (thẩm nhập gene) qua **backcrossing** (lai ngược liên tiếp). Đây là hiện tượng **tự nhiên và phổ biến** giữa các loài gần nhau; đó là cách đọc mặc định và phải nêu trước mọi diễn giải khác. [CONFIRMED — R11 mục 2]

Hai đại lượng suy giảm theo hai cách khác nhau, và chính chỗ này là mấu chốt:

- **Tỉ lệ tổ tiên** từ nguồn cho.
- **Chiều dài đoạn tổ tiên** liên tục trên nhiễm sắc thể, bị tái tổ hợp cắt nhỏ dần qua mỗi đời.

Đại lượng thứ hai là "đồng hồ": đoạn dài nghĩa là sự kiện lai gần đây, đoạn ngắn nghĩa là xa xưa. [CONFIRMED — R11, dẫn Gravel 2012 và Corbett-Detig & Nielsen 2017]

## 2. Phát hiện quan trọng nhất: hai kịch bản, hai số phận trái ngược

Mô phỏng cho thấy **cùng một sự kiện lai ban đầu 50/50 dẫn tới hai kết cục hoàn toàn khác nhau**, tùy quần thể lai sau đó làm gì. [CONFIRMED — `pipeline/p05_admixture/example/`, 37 test đạt, có test riêng phân biệt hai chế độ]

**Kịch bản A, lai ngược liên tiếp** (đúng hình ảnh "pha thêm nước" của chủ dự án). Mỗi đời cá thể lai ghép với một cá thể bản địa thuần:

| Đời | Tỉ lệ tổ tiên gốc | Số đoạn |
|---|---|---|
| 0 | 0,500 | 480 |
| 2 | 0,120 | 232 |
| 4 | 0,024 | 84 |
| 6 | 0,0069 | 35 |
| 10 | 0,00054 | 4 |

Phần gốc **biến mất**, xuống dưới 1% ngay ở đời thứ sáu.

**Kịch bản B, quần thể lai tự giao phối trong nội bộ** (tức nó trở thành một quần thể riêng, không còn lai ngược ra ngoài):

| Đời | Tỉ lệ tổ tiên gốc | Số đoạn | Trung vị chiều dài đoạn | Đoạn dài nhất |
|---|---|---|---|---|
| 0 | 0,500 | 1.200 | 100 cM | 100 cM |
| 10 | 0,474 | 4.293 | 18,97 cM | 100 cM |
| 40 | 0,491 | 11.330 | 6,73 cM | 100 cM |
| 80 | 0,506 | 18.580 | 3,86 cM | 75,9 cM |
| 200 | 0,496 | 28.085 | 2,17 cM | 62,0 cM |

Phần gốc **không hề biến mất**. Sau 200 đời nó vẫn chiếm gần đúng một nửa bộ gene. Cái thay đổi là cấu trúc: từ 1.200 đoạn dài thành 28.085 mảnh vụn.

**Đây là điều chỉnh quan trọng cho hình ảnh nước đường.** Đường chỉ loãng đi nếu ta liên tục thêm nước. Nếu pha một lần rồi đậy nắp thì đường vẫn còn nguyên trong cốc, chỉ tan đều ra. Sinh học cũng vậy: một quần thể lai trở thành loài riêng thì phần gốc ở lại vĩnh viễn, chỉ bị xé nhỏ.

## 3. Hệ quả cho cách đi tìm

Nếu giả thuyết là "chim hiện đại mang một bộ gene gốc từ nguồn khác, và quần thể đó về sau tự sinh sản", thì **cách tìm sai là đi tìm một tỉ lệ nhỏ lạ lẫm**. Phần gốc có thể chiếm tỉ lệ lớn và vẫn hoàn toàn vô hình với cách tìm đó, vì nó không nằm ở một chỗ mà nằm rải khắp nơi.

Cách tìm đúng là nhìn vào **cấu trúc đoạn**: có hay không những khối liên tục dài bất thường, mang tổ hợp biến thể khác hẳn phần còn lại của bộ gene. Đó chính là nguyên lý của nhóm phương pháp phát hiện thẩm nhập từ dòng dõi chưa có mẫu. [STRONG_INFERENCE — suy ra từ mục 2 cộng nguyên lý phương pháp ở R11]

Đáng chú ý: đoạn dài nhất vẫn còn 62 cM sau 200 đời. Vài đoạn dài luôn sống sót lâu hơn trung vị rất nhiều, và chính chúng là thứ các phương pháp bắt được.

## 4. Trạng thái thực tế của phương pháp

- Nhóm phương pháp **không cần bộ gene của nguồn cho** đã tồn tại và có một làn sóng công cụ mới trong năm nay. [CONFIRMED — R11]
- Một điểm hay bị hiểu lầm: không phải mọi công cụ "ghost" đều thật sự không cần nguồn; ít nhất một công cụ phổ biến vẫn cần một bộ gene cổ thật làm đầu vào. [CONFIRMED — R11]
- **Thẩm nhập từ dòng dõi ma đã được xác nhận ở chim**, ít nhất hai ca, cả hai đều thuộc bộ Sẻ như chào mào; một ca đã bình duyệt. [CONFIRMED — R11]
- **Chưa có nghiên cứu nào thuộc loại này cho họ Chào mào hay chi *Pycnonotus*.** Đây là khoảng trống thật. [CONFIRMED — R11 mục 8]

## 5. Cái chưa có, và vì sao nó chặn bước tiếp theo

**Không tìm được con số ngưỡng phát hiện đã qua kiểm chứng**: tỉ lệ tối thiểu bao nhiêu, đoạn ngắn tới đâu, sự kiện xa bao nhiêu đời thì các phương pháp mất khả năng. R11 gắn nhãn UNVERIFIED cho mọi con số cụ thể. [CONFIRMED — R11 mục 5 của TL;DR]

Hệ quả trực tiếp: lệnh `detect` trong mô phỏng nhận ngưỡng do người dùng đặt, nên **mọi kết luận từ nó hiện chỉ mang tính minh họa**, không phải phát biểu về khả năng phát hiện thật. Muốn có số thật phải đọc tài liệu phương pháp của từng công cụ, hoặc tự chạy hiệu chuẩn.

Giới hạn của chính mô phỏng: bỏ qua chọn lọc, bỏ qua thay đổi kích thước quần thể, bỏ qua biến thiên tốc độ tái tổ hợp dọc nhiễm sắc thể, giả định bản đồ di truyền đồng nhất. Nó cho **bậc độ lớn**, không cho con số dùng để công bố. Kịch bản B còn giả định toàn bộ quần thể khởi đầu ở tỉ lệ 0,5, tức một bầy lai thuần túy không còn cá thể bản địa nào.

## 6. Việc kế tiếp nếu đi hướng này

1. Lấy ngưỡng phát hiện thật từ tài liệu của từng công cụ, rồi chạy lại `detect` để có phát biểu có nghĩa.
2. Kiểm tra dữ liệu tối thiểu cần có cho *P. jocosus*: số cá thể, độ sâu giải trình tự. Bộ gene hiện tại chỉ ở mức scaffold, gần như chắc chắn không đủ.
3. Trước khi bàn tới nguồn lạ, phải loại trước khả năng thông thường và phổ biến hơn nhiều: thẩm nhập **giữa các loài chim gần nhau**, ví dụ giữa các loài cùng chi *Pycnonotus* vốn phân bố chồng lấn ở Việt Nam.
