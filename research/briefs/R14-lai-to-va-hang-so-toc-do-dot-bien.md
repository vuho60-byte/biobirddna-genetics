# R14 — Lại tổ (atavism) ở chim và các hằng số cần cho phép tính xác suất

Đọc `research/briefs/_COMMON.md`. Output: `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md`.

## Câu hỏi thật đằng sau brief này

Chủ dự án hỏi: **xác suất một gene cổ hiện lại là bao nhiêu phần trăm, và có bằng xác suất xuất hiện một đột biến mới không?**

Câu này tính được, nhưng chỉ khi có hằng số thật. Việc của bạn là mang về những hằng số đó kèm nguồn, **không** tự tính hộ.

Điểm mấu chốt phải làm rõ trong báo cáo: xác suất "hiện lại" phụ thuộc hoàn toàn vào **cách đặc điểm đó bị mất**. Ba trường hợp có xác suất chênh nhau hàng chục bậc độ lớn:

- (A) chương trình phát triển còn nguyên, chỉ một công tắc điều hòa bị tắt;
- (B) gene bị bất hoạt bởi nhiều đột biến tích lũy (giả gene);
- (C) gene bị xóa hẳn khỏi bộ gene.

## Câu hỏi phải trả lời

1. **Tốc độ đột biến ở chim**: số đột biến mới trên mỗi vị trí trên mỗi thế hệ, đo trực tiếp bằng giải trình tự bộ ba cha-mẹ-con. Ưu tiên các công trình đo thật ở chim (ví dụ dòng công trình trên *Ficedula*, *Taeniopygia*, gà). Cho **con số kèm khoảng tin cậy và loài đo được**. Nếu có nhiều ước lượng khác nhau, liệt kê hết, đừng chọn một.

2. **Tỉ lệ chuyển đổi so với đảo chuyển** (transition/transversion ratio) ở chim. Cần vì xác suất một base đổi thành **đúng một base cụ thể** không bằng một phần ba tốc độ chung.

3. **Thời gian thế hệ** của chim sẻ nhỏ trong tự nhiên, và của chào mào nếu có số riêng.

4. **Kích thước quần thể hiệu dụng** (Ne) điển hình của chim sẻ hoang dã, và của các loài nuôi nhốt nếu có số.

5. **Ca lại tổ đã ghi nhận ở chim**, với cơ chế phân tử nếu đã biết. Bắt buộc có: đột biến `talpid2` ở gà làm phôi mọc răng (Harris và cộng sự, khoảng 2006) — cơ chế chính xác là gì, răng đó có men không, phôi có sống được không. Các ca khác: ngón chân thừa, lông ở chân, bất kỳ ca nào khác trong tài liệu bình duyệt.

6. **Số đột biến bất hoạt trong gene răng của chim**: ENAM, AMBN, AMELX, AMTN, DSPP, MMP20. Với mỗi gene, ở các loài chim đã khảo sát: có bao nhiêu đột biến làm hỏng (dịch khung, mã dừng sớm, mất đoạn), và chúng có **chung nguồn gốc** ở các loài không. Đây là con số quyết định cho trường hợp (B). Nguồn chính có thể là Meredith và cộng sự 2014 *Science* (R12 đã dùng) — hãy vào bảng phụ lục lấy số cụ thể nếu mở được.

7. **Định luật Dollo**: phát biểu gốc, và các công trình hiện đại kiểm định nó bằng phát sinh chủng loại. Có ca phản ví dụ nào được chấp nhận (đặc điểm phức tạp mất rồi lấy lại được) không? Ví dụ hay được nhắc: ngón chân ếch, cánh côn trùng que.

8. **Ưu thế lai và sự suy giảm của nó**: ưu thế lai đời F1 có thật bao nhiêu, và hiện tượng suy giảm ở F2 (outbreeding depression / hybrid breakdown) được ghi nhận thế nào. Cần vì chủ dự án giả định quần thể lai giữ được sức khỏe vượt trội qua nhiều đời.

9. **Thẩm nhập thích nghi**: ca đã được chấp nhận rộng rãi nhất là gì. Ví dụ hay được nhắc là allele EPAS1 thích nghi độ cao ở người Tây Tạng đến từ người Denisova — xác minh và nêu chính xác. Có ca tương đương ở chim không?

## Ranh giới bắt buộc
- **Không tự tính xác suất.** Chỉ mang hằng số và nguồn về; phép tính do Main làm.
- Mỗi con số phải kèm loài đo được và phương pháp đo. Số ngoại suy từ loài khác phải ghi rõ.
- Nhãn C/S/U. Perplexity tối đa 3 lượt.
