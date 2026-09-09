# R11 — Mô hình "pha loãng": phát hiện dấu vết lai từ một nguồn chưa biết

Đọc `research/briefs/_COMMON.md` (luật nhãn, format 6 mục, giới hạn Perplexity). Output: `research/raw/R11-introgression-ghost-lineage.md`.

## Câu hỏi thật đằng sau brief này

Chủ dự án mô tả một mô hình bằng hình ảnh pha nước đường: một bộ gene "gốc" trộn với một bộ gene "bản địa", rồi qua nhiều đời lai ngược thì phần gốc loãng dần cho tới khi hòa vào nền. Anh ấy hỏi liệu có thể tái dựng lại quá trình đó không.

Mô hình này **trùng khít** với một khái niệm chuẩn của di truyền quần thể: **introgression** (thẩm nhập gene) qua **backcrossing** (lai ngược liên tiếp). Và có một nhánh phương pháp trả lời đúng câu hỏi khó nhất của anh ấy: phát hiện phần đóng góp của một dòng dõi **không có mẫu, không có bộ gene tham chiếu** — gọi là **ghost introgression** hoặc **archaic introgression**.

Việc của bạn là thu thập những gì đã biết về nhóm phương pháp này, và đánh giá xem có áp được cho chim hay không.

## Câu hỏi phải trả lời

1. **Toán của sự loãng dần.** Tỉ lệ tổ tiên từ nguồn cho sau F1, F2, ... Fn dưới lai ngược là bao nhiêu (công thức)? Quan trọng hơn: **phân bố chiều dài đoạn tổ tiên** (ancestry tract length) suy giảm thế nào theo số đời do tái tổ hợp? Nêu công thức và nguồn giáo trình hoặc bài báo gốc. Đây là "đồng hồ" cho phép ước lượng sự kiện lai xảy ra cách đây bao nhiêu đời.

2. **Phương pháp phát hiện khi CÓ bộ gene nguồn**: thống kê D (ABBA-BABA), f4-ratio, f_dM, RND, QuIBL. Mỗi cái đo gì, giả định gì, hỏng khi nào.

3. **Phương pháp phát hiện khi KHÔNG CÓ bộ gene nguồn** (phần quan trọng nhất của brief): S*, Sprime, ArchIE, IBDmix, và bất kỳ phương pháp nào khác bạn tìm được. Với mỗi cái: nguyên lý, dữ liệu tối thiểu cần có, đã dùng cho loài nào, tỉ lệ dương tính giả đã báo cáo.

4. **Giới hạn phát hiện.** Với tỉ lệ đóng góp thấp tới đâu, và sự kiện lai xa tới bao nhiêu đời, thì các phương pháp trên còn phát hiện được? Có bài báo nào đo trực tiếp ngưỡng này bằng mô phỏng không? Đây là câu quyết định: nó cho biết một sự kiện lai cổ xưa **có thể trở nên vô hình** hay không.

5. **Introgression đã ghi nhận ở CHIM**: loài nào, tỉ lệ bao nhiêu, phương pháp gì. Ưu tiên chim sẻ (Passeriformes). Có ca nào ở Pycnonotidae hoặc chi *Pycnonotus* không?

6. **Có công trình nào công bố phát hiện đóng góp từ một dòng dõi ma (ghost/unsampled) trong bộ gene chim không?** Nếu có, nêu chính xác; nếu không tìm thấy, ghi rõ "không tìm thấy" — đó là câu trả lời hợp lệ và có giá trị.

7. **Yêu cầu dữ liệu thực tế.** Để chạy một phân tích ghost introgression cho *P. jocosus* cần gì: bao nhiêu cá thể, độ sâu giải trình tự bao nhiêu, có cần bộ gene mức nhiễm sắc thể không, có cần bản đồ tái tổ hợp không. Ước lượng chi phí nếu tìm được nguồn.

8. **Dữ liệu tải được**: phần mềm (URL, license), bộ dữ liệu chim đa cá thể đã công bố có thể tái dùng.

## Ranh giới bắt buộc

- Trình bày trung lập. Introgression là hiện tượng **tự nhiên phổ biến** (lai giữa các loài gần nhau); đó là cách đọc mặc định và phải nêu trước.
- **Không** kết luận gì về nguồn gốc ngoài Trái Đất. Brief này thu thập phương pháp; việc diễn giải thuộc về Main và phải dựa trên số liệu của chính dự án.
- Nhãn C/S/U như thường lệ. Perplexity tối đa 3 lượt.
