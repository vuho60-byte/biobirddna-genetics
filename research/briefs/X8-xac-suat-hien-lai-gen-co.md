# X8 — Mô hình phân tử: xác suất một gene cổ hiện lại

Chỉ tạo/sửa trong `pipeline/p06_atavism/**` và MỘT dòng thêm vào `tasks/primordial-dna-handoff.md`. **Không có lệnh xóa nào** (rm/rmdir/del/Remove-Item) trong bất kỳ script nào. Python 3.11 **stdlib**. Không numpy/scipy, không docker, không tải dữ liệu.

## Câu hỏi khoa học

Chủ dự án hỏi hai câu, và chúng có câu trả lời bằng số:

1. Xác suất một gene cổ **hiện lại** là bao nhiêu?
2. Có bằng xác suất **xuất hiện một đột biến mới** không?

Chủ dự án đề nghị mô hình hóa ở mức phân tử: một locus gồm các phần tử đặt tên **A, B, C, ... Z**; phần tử **A** là phần tử đã bị đột biến làm hỏng chức năng, các phần tử còn lại còn nguyên. Đây chính là cách đóng khung đúng, vì câu trả lời phụ thuộc hoàn toàn vào **có bao nhiêu phần tử phải quay ngược lại**.

## Mô hình phải cài

Một locus = danh sách phần tử có tên A, B, C, ... Z. Mỗi phần tử có trạng thái `goc` hoặc `dot_bien`. Ban đầu k phần tử đầu tiên (A, hoặc A và B, hoặc A B C...) ở trạng thái `dot_bien`; phần còn lại `goc`.

**Xác suất quay ngược một phần tử cụ thể trong một thế hệ, trên một bản sao nhiễm sắc thể:**

    p_rev = mu * p_specific

trong đó `mu` là tốc độ đột biến trên mỗi vị trí trên mỗi thế hệ, và `p_specific` là xác suất đổi thành **đúng** base gốc chứ không phải hai base kia. Với tỉ lệ tốc độ chuyển đổi trên đảo chuyển là `kappa`, mỗi base có 1 đường chuyển đổi và 2 đường đảo chuyển:

    p_specific = kappa / (kappa + 2)      nếu base gốc nằm ở đường chuyển đổi
    p_specific = 1 / (kappa + 2)          nếu ở một đường đảo chuyển cụ thể

Cài cả hai, để người dùng chọn bằng cờ; mặc định lấy trường hợp thuận lợi nhất (chuyển đổi) để **con số đưa ra là giới hạn TRÊN của khả năng hiện lại**. Ghi rõ lựa chọn này trong README.

**Hai con đường quay ngược, phải cài cả hai:**

- **Đồng thời**: cả k phần tử cùng quay ngược trong một giao tử. Xác suất mỗi giao tử mỗi thế hệ = tích các `p_rev`. Số ca kỳ vọng mỗi thế hệ trong quần thể = `2 * N * tích(p_rev)`.
- **Tuần tự**: các phần tử quay ngược lần lượt, mỗi bước trung tính (gene đã hỏng sẵn nên phần tử quay ngược một mình không cho lợi thế nào) và phải cố định trong quần thể trước khi bước sau xảy ra. Đột biến trung tính cố định với xác suất `1/(2N)`, nên thời gian chờ mỗi bước ≈ `1 / (mu * p_specific)` thế hệ, độc lập với N. Tổng ≈ `k / (mu * p_specific)`. Ghi rõ đây là kết quả chuẩn của di truyền quần thể trung tính và **nếu bạn không xác minh được thì gắn nhãn giả định trong README**, đừng khẳng định.

Đường tuần tự cho thời gian NGẮN hơn đường đồng thời khi k ≥ 2, nên nó là kịch bản thuận lợi nhất và phải là con số chính đưa ra.

**So sánh bắt buộc — đây là câu hỏi 2 của chủ dự án.** Xác suất xuất hiện **một đột biến mới bất kỳ** ở đâu đó trong locus dài `L` base:

    p_new = 1 - (1 - mu)^L  ≈ L * mu

Tỉ số `p_new / p_rev` cho biết đột biến mới dễ xảy ra hơn bao nhiêu lần so với quay ngược đúng một phần tử cụ thể. Đây là con số trả lời thẳng câu hỏi.

## Việc — `pipeline/p06_atavism/reversion_prob.py`, bốn lệnh con

1. `locus --elements A,B,C --mutated A --mu 4.6e-9 --kappa 2.0 --n-diploid 1e6 --generation-years 2 [--transversion]`
   In bảng trạng thái từng phần tử, rồi xác suất và thời gian chờ theo cả hai con đường, tính bằng thế hệ VÀ bằng năm.

2. `compare --locus-length 1000 --mu ... --kappa ... --mutated-count 1`
   So `p_new` với `p_rev`, in tỉ số. Đây là lệnh trả lời câu hỏi 2.

3. `sweep --mu-list ... --k-list 1,2,3,5 --n-list ... --out sweep.tsv`
   Bảng thời gian chờ theo lưới tham số. Cần vì hằng số thật còn đang được brief R14 tra; kết quả phải bền vững trước khoảng giá trị chứ không phụ thuộc một con số.

4. `report --sweep sweep.tsv --out report.md`
   Bảng gọn kèm phần **cảnh báo bắt buộc**: (a) "xuất hiện" khác "nhìn thấy được" khác "lan ra quần thể" — một bản quay ngược ở thể dị hợp sẽ không biểu hiện nếu tính trạng là lặn; (b) mô hình bỏ qua chọn lọc, bỏ qua đột biến mất đoạn (không quay ngược được bằng con đường nào), bỏ qua chuyển đổi gene; (c) nếu phần tử bị **xóa hẳn** thì xác suất quay ngược bằng 0 theo con đường này, không phải một số nhỏ.

## Test (`pipeline/p06_atavism/tests/test_reversion_prob.py`, unittest)

- `p_specific` với hai lựa chọn chuyển đổi/đảo chuyển cộng lại đúng theo tỉ lệ kappa; kappa = 1 cho đúng 1/3 cho mọi đường.
- Thời gian chờ đường tuần tự tỉ lệ tuyến tính với k và **không** phụ thuộc N; đường đồng thời thì phụ thuộc N. Đây là test phân biệt hai con đường, bắt buộc phải có.
- Tăng `mu` làm giảm thời gian chờ đúng tỉ lệ nghịch.
- `compare` với locus 1 base và 1 phần tử đột biến cho tỉ số đúng bằng `1/p_specific`.
- k = 0 (không phần tử nào hỏng) phải báo lỗi rõ ràng, không trả về vô cực im lặng.
- Số rất nhỏ không bị tràn xuống 0: dùng `math.log` khi cần, kiểm bằng k = 10.

## Tiêu chí xong
- `python -B -m unittest discover -s pipeline/p06_atavism/tests -v` sạch; `--help` mỗi lệnh con thoát 0.
- Chạy thật `sweep` với mu trong khoảng 1e-9 đến 1e-8, k trong 1,2,3,5, N trong 1e4,1e6,1e8; ghi vào `pipeline/p06_atavism/example/`.
- README nêu công thức, ý nghĩa cột, ba cảnh báo trên, và ghi rõ hằng số mặc định là **tạm**, chờ `research/raw/R14-*.md`.
- Thêm đúng 1 dòng `tasks/primordial-dna-handoff.md`.

## Trả về ≤10 dòng
file tạo · dòng tổng kết unittest · thời gian chờ (năm) cho k=1 và k=3 ở mu=4,6e-9 · tỉ số p_new/p_rev với locus 1000 base · xác nhận không lệnh xóa · giả định chính.
