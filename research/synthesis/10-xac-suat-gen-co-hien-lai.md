# 10 — Xác suất một gene cổ hiện lại: con số

Ngày: 2026-09-09. Mô hình: `pipeline/p06_atavism/reversion_prob.py` (41/41 test đạt, tôi tự chạy lại và tự kiểm phép tính bằng tay). Hằng số: `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md`.

## 1. Hằng số dùng trong phép tính

| Hằng số | Giá trị | Nguồn | Mức xác minh |
|---|---|---|---|
| Tốc độ đột biến mỗi vị trí mỗi thế hệ | 4,6 × 10⁻⁹ | *Ficedula albicollis*, Smeds, Qvarnström & Ellegren 2016, *Genome Research*, [DOI](https://doi.org/10.1101/gr.204669.116) | CONFIRMED — Main tự mở PubMed xác nhận |
| Tỉ lệ chuyển đổi trên đảo chuyển | 2,67 (KTC 95%: 1,56–5,44) | cùng bài, 47 chuyển đổi trên 8 đảo chuyển trong 55 đột biến mới | CONFIRMED — R14 đọc toàn văn PMC |
| Số năm mỗi thế hệ | 2 | giả định làm việc cho chim sẻ nhỏ | GIẢ ĐỊNH — chưa có số riêng cho chào mào |
| Kích thước quần thể lưỡng bội | 10⁶ | giả định làm việc | GIẢ ĐỊNH |

Chọn đường **chuyển đổi** (dễ xảy ra hơn đảo chuyển) nên mọi con số dưới đây là **giới hạn TRÊN của khả năng hiện lại**. Thực tế khó hơn.

## 2. Câu hỏi 1 — xác suất hiện lại

Phải tách hai câu hỏi khác nhau mà tiếng Việt thường gộp làm một:

- **Xuất hiện**: có ít nhất một cá thể mang bản quay ngược.
- **Trở thành đặc điểm của loài**: bản quay ngược lan ra và cố định trong quần thể.

Hai con số cách nhau sáu bậc độ lớn.

### Trường hợp k = 1 (một phần tử hỏng, phần còn lại nguyên vẹn)

| Đại lượng | Giá trị |
|---|---|
| Xác suất mỗi bản sao mỗi thế hệ | 2,63 × 10⁻⁹ |
| Số ca kỳ vọng mỗi thế hệ trong quần thể 10⁶ | 0,0053 |
| Chờ tới khi **xuất hiện** ở một cá thể | 190 thế hệ ≈ **380 năm** |
| Chờ tới khi **cố định** trong loài | 3,80 × 10⁸ thế hệ ≈ **760 triệu năm** |

### Trường hợp k = 3 (ba phần tử cùng hỏng)

| Đại lượng | Giá trị |
|---|---|
| Xác suất mỗi giao tử mỗi thế hệ (cùng lúc) | 1,82 × 10⁻²⁶ |
| Chờ theo đường cùng lúc | 5,50 × 10¹⁹ năm |
| Chờ theo đường tuần tự (nhanh hơn) | **2,28 × 10⁹ năm** |

Để so sánh: chim hiện đại tồn tại khoảng 10⁸ năm; vũ trụ khoảng 1,4 × 10¹⁰ năm.

**Trả lời bằng phần trăm như chủ dự án hỏi:** với một phần tử, khoảng 0,0000005% mỗi cá thể mỗi thế hệ. Với ba phần tử, con số nhỏ tới mức viết ra không còn nghĩa.

## 3. Câu hỏi 2 — có bằng xác suất đột biến mới không?

**Không.** Với locus dài 1.000 base và một phần tử cần quay ngược:

| Đại lượng | Giá trị |
|---|---|
| Xác suất xuất hiện một đột biến MỚI bất kỳ trong locus | 4,60 × 10⁻⁶ |
| Xác suất quay ngược đúng phần tử đó | 2,63 × 10⁻⁹ |
| **Tỉ số** | **1.749 lần** |

Đột biến mới dễ xảy ra hơn khoảng 1.750 lần, vì nó được rơi vào bất kỳ chỗ nào trong nghìn vị trí, còn quay ngược phải trúng đúng một chỗ và đúng một base. Tỉ số này tăng theo lũy thừa với số phần tử phải quay ngược.

## 4. Vì sao gà talpid2 mọc răng nhưng răng không có men

Đây là phép thử của chính tự nhiên, và nó khớp đúng khung ba trường hợp.

- Đột biến talpid2 nằm ở gene C2CD3, mất 19 base gây mã dừng sớm, làm rối loạn tín hiệu Hedgehog và đổi vị trí trung tâm tín hiệu phát triển. **Không mất mã di truyền nào cho răng.** Đây là trường hợp k nhỏ. [CONFIRMED — Chang và cộng sự 2014, *Development*, DOI 10.1242/dev.105924, qua R14]
- Nhưng sáu gene cấu trúc răng (ENAM, AMBN, AMELX, AMTN, DSPP, MMP20) mang đột biến bất hoạt **dùng chung** ở gần như mọi bộ chim, khớp một lần mất duy nhất khoảng 116 triệu năm trước. Đây là trường hợp k lớn. [CONFIRMED — Meredith và cộng sự 2014, *Science*, qua R12 và R14]

Bật lại được công tắc vị trí, không sửa được gene đã hỏng. Kết quả: có răng, không có men. [STRONG_INFERENCE — kết hợp hai nguồn CONFIRMED, không phải trích trực tiếp]

Và điều này giải thích luôn nghịch lý biểu kiến: talpid2 tồn tại như một đột biến hiếm trong phòng thí nghiệm (thời gian chờ **xuất hiện** ngắn) nhưng chim không có răng (thời gian chờ **cố định** dài hơn tuổi của cả lớp Chim).

## 5. Đính chính quan trọng: định luật Dollo không tuyệt đối

Hai bài tổng quan độc lập 2024–2025 đều kết luận định luật Dollo **không tuyệt đối**. Nó đúng khi không còn cấu trúc tương đồng nào sót lại, nhưng có thể bị phá khi cấu trúc tương đồng nối tiếp vẫn tồn tại ở nơi khác trên cơ thể. [CONFIRMED — Elmer & Clobert 2024; Jackman, Lynch & Gibert 2025, qua R14]

Ca phản ví dụ được cả hai dẫn: ếch *Gastrotheca guentheri* mất răng hàm dưới ít nhất 230 triệu năm rồi **lấy lại** khoảng 5–17 triệu năm trước, vì nó vẫn giữ răng hàm trên nên bộ máy làm răng chưa bao giờ ngừng chịu chọn lọc. [CONFIRMED — Wiens 2011, qua R14]

**Hệ quả cho dự án:** thứ quyết định khả năng quay lại không phải thời gian đã trôi qua, mà là **bộ máy có còn được duy trì ở đâu đó không**. Ở chim, bộ máy làm men răng không được duy trì ở bất kỳ đâu. Nên nó đi hẳn.

Lưu ý: ca cánh côn trùng que (Whiting 2003) hay được nhắc như phản ví dụ thì bị chính Wiens 2011 ghi nhận là có thể do lỗi phương pháp. Không dùng làm bằng chứng.

## 6. Giới hạn của phép tính

- Thời gian chờ đường tuần tự dùng kết quả chuẩn của lý thuyết trung tính. Đây là **giả định chưa tự mở nguồn xác minh** trong dự án, ghi rõ trong README của mô hình.
- Mô hình bỏ qua chọn lọc, bỏ qua chuyển đổi gene, bỏ qua đột biến mất đoạn.
- Nếu phần tử bị **xóa hẳn** thì xác suất quay ngược bằng 0 theo con đường này, không phải một số nhỏ.
- Một bản quay ngược ở thể dị hợp sẽ không biểu hiện nếu tính trạng là lặn, nên "xuất hiện" vẫn chưa bằng "nhìn thấy được".
- Chưa có số đo tốc độ đột biến trực tiếp cho gà, và chưa có thời gian thế hệ riêng cho chào mào.
