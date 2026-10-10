# 05 — Mô tả phát hiện và hướng mở rộng

Ngày 2026-09-09. Mọi số trong tài liệu này tái lập được từ `data/a2/`, `data/a2_accel/`, `data/anomaly/`.

---

## PHẦN A — Mô tả kỹ ba kết quả đã có

### A1. Đối chứng âm: quy trình phân biệt được vùng chịu ràng buộc với vùng tự do biến đổi

**Thiết kế.** Chạy **cùng một quy trình, cùng tham số, cùng ngày** trên hai bộ dữ liệu đối lập, đều lấy từ phép đo phyloP trên alignment 363 loài chim:
- bộ **conserved**: vùng thay thế base **chậm hơn** mô hình trung tính (FDR < 5%);
- bộ **accelerated**: vùng thay thế base **nhanh hơn** mô hình trung tính (cùng FDR).

Hai bộ khác nhau đúng một điều: dấu của độ lệch so với trung tính. Mọi thứ khác giống hệt.

**Kết quả.**

| | Bảo tồn | Tiến hóa nhanh |
|---|---|---|
| Phần tử đưa vào | 2.947.588 | 4.054.556 |
| Base gốc trên bộ gene gà | 134.568.082 (12,85%) | 89.349.479 (8,53%) |
| Ánh xạ được sang chào mào | **248.461** | **1.133** |
| Tỉ lệ ánh xạ | **8,43%** | **0,03%** |

Chênh **281 lần**, dù bộ tiến hóa nhanh có **nhiều phần tử hơn** (4,05 so với 2,95 triệu).

**Ý nghĩa — nói cho đúng.** Đây **không phải phát hiện sinh học mới**. Đây là **kiểm chứng phương pháp**, và nó đóng lại nghi ngờ lớn nhất về toàn bộ kết quả: rằng 6% bộ gene tìm được chỉ là thứ công cụ căn trình tự tình cờ khớp. Nếu vậy, bộ tiến hóa nhanh cũng phải khớp tương đương. Nó không.

Kết quả này đồng thời là một **thang đo** cho mọi kết quả về sau: bất kỳ tập vùng nào ánh xạ được ở mức gần 0,03% là nhiễu; gần 8% là tín hiệu thật.

**Giới hạn.** Vì gần như không gene nào có vùng tiến hóa nhanh ánh xạ được, **không** tính được tỉ số bảo tồn/tăng tốc theo từng gene như dự tính ban đầu. Đối chứng âm này có giá trị ở cấp **toàn quy trình**, không ở cấp gene. Muốn xếp hạng gene vẫn phải dùng xáo trộn vị trí ngẫu nhiên (chưa làm).

### A2. Bản đồ lõi bất biến đầu tiên cho họ Chào mào

- 248.461 vùng, **61,50 Mb sau khi gộp trùng lặp = 6,00% bộ gene chào mào**.
- **59,8%** toàn bộ vùng mã hóa protein của chào mào nằm trong lõi; làm giàu **9,96 lần** so với rải ngẫu nhiên (giao base-level).
- Phân bố lõi: 19,2% trong vùng mã hóa, 33,3% trong thân gene (gồm intron), phần còn lại ở vùng giữa các gene.
- Xếp gene theo **mật độ** lõi trên mỗi kb locus đưa **cụm Hoxa** (Hoxa3, Hoxa5, Hoxa6, Hoxa11) lên đầu, trên 1.800 base lõi mỗi kb. Xếp theo tổng base thì đứng đầu là gene dài (Znf521, Sox6, Foxp1) — hai bảng không giao nhau một gene nào.

**Tính mới:** bảo tồn sâu ở vùng mã hóa và ở cụm Hox là kiến thức chuẩn. Cái mới là **bộ dữ liệu cho chính họ Pycnonotidae**, chưa từng được công bố, tái lập được bằng mã trong repo.

### A3. Không có dấu vết đoạn DNA ngoại lai trong vùng mã hóa

Quét 14.003 vùng mã hóa bằng ba thước đo thành phần độc lập (GC/GC3, cách dùng codon, phổ 4-mer):

| Số cờ | Số trình tự |
|---|---|
| 3 cờ (nghi ngờ mạnh) | **0** |
| 2 cờ | 27 |
| 1 cờ | 89 |
| 0 cờ | 13.887 (99,17%) |

27 trình tự hai cờ đều ngắn (324–549 bp), tỉ lệ GC bình thường, chỉ lệch ở phổ 4-mer — chữ ký của **trình tự ngắn thống kê kém tin cậy hoặc vùng lặp**, không phải của DNA khác dòng dõi. Kết quả **âm tính**, và âm tính ở đây là kết quả có giá trị: nó thu hẹp không gian giả thuyết bằng số liệu.

---

## PHẦN B — Mười hướng mở rộng, xếp theo chi phí và giá trị

### Nhóm 1 — Hoàn thiện cái đang dở (rẻ, làm được ngay)

**B1. Căn đoạn ngắn bằng LAST.** Hiện nhóm dưới 100 bp chiếm ~78% số phần tử nhưng chỉ ánh xạ được 0,25–4,7%. Mã đã viết, 14 phép thử đạt, chưa chạy. Kỳ vọng nâng độ phủ từ 6% lên có thể 8–10% bộ gene và **sửa được thiên lệch** hiện tại (bức tranh đang nghiêng về vùng mã hóa vì đoạn dài dễ căn hơn).
*Chi phí: một buổi chạy. Giá trị: cao — nó quyết định con số phân bố có tin được không.*

**B2. Kiểm định nền bằng xáo trộn vị trí.** Viết bộ xáo trộn ngẫu nhiên bằng Python thuần: giữ nguyên số lượng và phân bố độ dài phần tử lõi, rải lại ngẫu nhiên trên bộ gene (tránh khoảng trống lắp ráp), lặp 1.000 lần, so số base lõi thật của từng gene với phân phối ngẫu nhiên. Đây mới là phép kiểm đúng cho câu "gene này giàu lõi hơn mong đợi".
*Chi phí: một ngày viết + vài giờ chạy. Giá trị: cao — biến bảng xếp hạng thành phát biểu có kiểm định.*

**B3. Tái dựng trình tự tổ tiên cho 20 locus.** Sản phẩm C đã lên kế hoạch, dữ liệu mất trong sự cố, cần chạy lại. Cho ra trình tự tổ tiên tại ba mốc và bảng vị trí khác biệt giữa tổ tiên và chào mào hiện đại.
*Chi phí: nửa ngày. Giá trị: cao — đây là "DNA gốc" theo nghĩa đen nhất mà tính toán cho được.*

### Nhóm 2 — Câu hỏi khoa học mới, khả thi với dữ liệu đang có

**B4. Vùng tăng tốc riêng của nhánh chào mào.** Dữ liệu UCSC đo tăng tốc trên **toàn cây chim**. Chúng ta có ba mô hình trung tính và có alignment chứa chào mào, nên chạy được phyloP ở chế độ so **một nhánh riêng** để tìm vùng bảo tồn ở mọi loài chim **nhưng thay đổi nhanh riêng ở nhánh chào mào**. Đây chính là phương pháp đã tạo ra khái niệm vùng tăng tốc ở người.
*Trả lời câu: **cái gì làm chào mào khác các loài chim khác?** Chưa ai làm cho họ này.*
*Chi phí: cần công cụ PHAST trong Docker và một cây loài; hai đến ba ngày. Giá trị: cao nhất trong danh sách về tính mới.*

**B5. Phép thử gene răng hỏng.** Vòng nghiên cứu trước đã xác nhận năm gene men răng không còn bản ghi nào ở gà. Quét chúng trong bộ gene chào mào bằng cách so protein cá sấu và rùa lên bộ gene. Nếu tìm thấy mảnh vỡ mang **cùng kiểu hỏng** như ở các bộ chim khác, đó là bằng chứng trực tiếp cho tổ tiên chung, kiểm được trên chính loài của bạn.
*Chi phí: cần công cụ tìm kiếm trình tự và bộ gene ngoài nhóm (đã có accession); hai ngày. Giá trị: cao — đây là phép thử phân biệt giả thuyết, không chỉ mô tả.*
*Trạng thái 2026-10-10: mã + test xong (`pipeline/p07_tooth_genes/`, brief X9); chưa chạy dữ liệu thật — cần duyệt image BLAST. Ghi chú: Meredith 2014 khảo sát 6 gene (năm gene trên + MMP20) theo R12 dòng #23; công cụ quét cả 6 + GULO, đối chứng dương MEPE/SPARC.*

**B6. So sánh lõi giữa các loài chim hót.** Chạy đúng quy trình này cho sẻ vằn và sẻ ngô lớn (bộ gene chất lượng cao hơn nhiều), rồi so ba bản đồ. Phần chung là lõi của chim hót; phần chỉ có ở chào mào là ứng viên cho đặc trưng riêng.
*Chi phí: mỗi loài một lần chạy, khoảng nửa ngày. Giá trị: trung bình cao; cũng là phép kiểm chéo cho chất lượng bộ gene chào mào vốn còn vụn.*

### Nhóm 3 — Nối thẳng vào việc nuôi chim

**B7. Lõi ở các gene sắc tố.** Kiểm xem MC1R, TYR, SLC45A2, OCA2, MLPH và vùng điều hòa quanh chúng được lõi phủ tới đâu. Gene bị phủ dày nghĩa là ít dung thứ đột biến, nên đột biến ở đó thường gây hậu quả nặng; gene phủ thưa thì dễ sinh biến thể sống được. Đây là cơ sở để đoán **loại đột biến màu nào khả dĩ xuất hiện** ở chào mào.
*Chi phí: vài giờ, dữ liệu đã có. Giá trị: cao cho ứng dụng thực tế, nối thẳng vào công cụ di truyền màu lông.*

**B8. Lõi ở gene liên quan tiếng hót.** Cùng cách làm với FoxP2 và các gene học hót. Trả lời câu hỏi thứ hai của bạn từ đầu dự án: phần "có sẵn trong DNA" của tiếng hót nằm ở đâu và bảo tồn tới mức nào.
*Chi phí: vài giờ. Giá trị: cao về mặt câu chuyện, và đo được.*

### Nhóm 4 — Dài hơi

**B9. Công bố bộ dữ liệu.** Bản đồ lõi đầu tiên cho họ Chào mào, kèm mã và mã băm, đủ chuẩn cho một bài dạng công bố dữ liệu. Cần hoàn thành B1, B2 và bổ sung đánh giá chất lượng bộ gene trước.

**B10. Bổ sung bộ gene chất lượng cao cho chào mào.** Hạn chế lớn nhất hiện nay là bộ gene chào mào ở mức scaffold, đoạn trung vị 218 kb. Mọi kết luận về trật tự gene và về vùng điều hòa xa đều bị chặn bởi điều này. Cách gỡ: chờ bộ gene mới, hoặc hợp tác với nhóm có mẫu.

---

## PHẦN C — Đề xuất thứ tự

1. **B1** (đoạn ngắn) — vì nó sửa thiên lệch của mọi con số hiện có.
2. **B7 + B8** (gene sắc tố và tiếng hót) — rẻ, trả lời đúng câu hỏi gốc của dự án, dùng được ngay.
3. **B2** (kiểm định nền) — biến bảng xếp hạng thành phát biểu có kiểm định.
4. **B5** (gene răng) — phép thử phân biệt giả thuyết đầu tiên thực sự chạy được.
5. **B4** (vùng tăng tốc riêng chào mào) — hướng có tính mới cao nhất, làm khi ba việc trên xong.

Ba việc đầu nằm gọn trong hai ngày làm việc với tài nguyên hiện tại.


---

## PHẦN D — Đính chính sau phản biện (2026-09-09, cùng ngày)

Ba phát biểu ở Phần A đã bị hạ mức sau khi Antigravity phản biện và tôi kiểm lại bằng chính tệp chú giải.

### D1. Bảng gene "giàu lõi" phần lớn là giả tạo do mô hình gene cụt

Kiểm tỉ lệ CDS trên chiều dài locus trong GFF chào mào:

| Gene | Locus (bp) | CDS (bp) | CDS/locus | Kết luận |
|---|---|---|---|---|
| Arc | 1.185 | 1.185 | 1,00 | mô hình gene CỤT — chỉ có CDS, không UTR/intron |
| Drd1 | 1.356 | 1.356 | 1,00 | mô hình gene CỤT |
| Ntrk2 | 1.116 | 327 | 0,29 | cụt nặng (gene thật ở chim dài hàng trăm kb) |
| Robo1 (1 trong 3 bản ghi) | 327 | 327 | 1,00 | mảnh vụn |
| Egr1 | 1.835 | 1.527 | 0,83 | ngắn, ít intron |
| Foxp2 | 91.384 | 1.818 | 0,02 | bình thường |
| Foxp1 | 160.237 | 2.112 | 0,01 | bình thường |

Vì phần mã hóa vốn giàu lõi nhất, gene có mẫu số bằng đúng CDS sẽ có mật độ thổi phồng. **Fold 17,7 của Drd1 và 17,0 của Arc gần như hoàn toàn là hiện tượng này.** Foxp2 (3,9) và Foxp1 (5,0) vẫn đáng tin vì mô hình gene đầy đủ.

**Sửa bắt buộc trước khi dùng bảng gene:** lọc bỏ gene có `CDS/locus > 0,9`, hoặc chuyển sang mô hình nền **phân tầng** (tính kỳ vọng riêng cho CDS, intron, vùng liên gene), hoặc dùng hoán vị đối sánh theo chiều dài CDS và số exon.

### D2. "Gene sắc tố ít bị ràng buộc" — không đứng vững

Sox10 (fold 4,7) và Kit (1,4) điều hòa mào thần kinh và tạo máu; đột biến mất chức năng thường gây chết phôi. Không thể gọi nhóm này là ít ràng buộc.

**Cách giải thích hợp lý hơn — thiên lệch quan sát:** người nuôi chim chỉ thấy đột biến ở những gene mà con vật còn sống để mang màu. Đột biến ở gene thiết yếu chết từ giai đoạn phôi nên không bao giờ vào danh mục. Điều này giải thích vì sao danh mục màu trong nghề tập trung vào một nhóm gene hẹp (OCA2, TYRP1, SLC45A2, MLPH) chứ không rải đều trên mọi gene sắc tố.

Ngoài ra Oca2 (fold 0,3) là locus lớn với intron rất dài, nên mật độ bị pha loãng cơ học — lại là lỗi mẫu số như D1, chiều ngược lại.

### D3. Đối chứng âm: giá trị hẹp hơn tôi đã viết

Vùng tăng tốc khó ánh xạ **vì chúng phân kỳ**, mà phân kỳ chính là định nghĩa của chúng. Một phần của chênh lệch 281 lần là hệ quả tất yếu của thuật toán căn trình tự, gần với lập luận vòng tròn.

- **Nó loại bỏ được:** giả thuyết "công cụ căn khớp bừa bất kỳ đoạn nào".
- **Nó KHÔNG chứng minh:** vùng bảo tồn có chức năng sinh học.
- **Cảnh báo thêm:** vùng tăng tốc chịu chọn lọc dương hoặc chuyển đổi gene thiên vị GC, nên không phải "chuỗi trung tính ngẫu nhiên" — nó không phải đối chứng âm chuẩn theo nghĩa thống kê.
- **Chưa kết luận được:** 91,57% vùng bảo tồn không ánh xạ được có thể do khoảng trống lắp ráp của bộ gene chào mào (mức scaffold), không nhất thiết do chào mào mất chúng.

**Đối chứng âm đúng nghĩa vẫn phải là hoán vị vị trí ngẫu nhiên (B2), chưa làm.**
