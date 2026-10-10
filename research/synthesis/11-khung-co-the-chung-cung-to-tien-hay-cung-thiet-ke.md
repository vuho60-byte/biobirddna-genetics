# 11 — Khung cơ thể chung: cùng tổ tiên hay cùng bản thiết kế? (tổng hợp R15)

Ngày: 2026-10-10. Nguồn chính: `research/raw/R15-khung-co-the-chung-thu-bac-long-nhau.md` (41 dòng phát hiện; theo `research/validate_raw.py`: 27 CONFIRMED / 8 STRONG_INFERENCE / 3 UNVERIFIED tính theo chuỗi khớp trong dòng, 4 dòng là ghi chú kỹ thuật không mang nhãn). Brief: `research/briefs/R15-khung-co-the-chung-thu-bac-long-nhau.md`. Đối chiếu chéo: R12 (gene răng chim), R13 (dấu hiệu bộ gene được thiết kế), R10 và `07-tra-loi-gia-thuyet-gene-la.md` (ERV ở chim), `docs/05` (kết quả A2 và quét thành phần của dự án).

Nhãn giữ nguyên của nguồn, không nâng. Bản này **không mở lại nguồn nào**; mọi nhãn là kế thừa từ R15. Số `#n` là số dòng trong bảng phát hiện của R15.

---

## 0. Câu hỏi của chủ dự án, nói lại cho đúng

Quan sát: gần như mọi động vật có xương sống dùng chung một bộ khung — đối xứng hai bên, cột sống, tim, gan, phổi; ở chim thêm túi khí, xương rỗng, cùng kiểu xương sườn. Suy ra: chúng được dựng từ một bộ khung có sẵn rồi sửa theo mục đích, nên "cấu trúc giống nhau thì cùng một DNA".

Quan sát **đúng**, và R15 lẫn `docs/00` Q1 đều ghi nhận như quan sát thật. Chỗ phải kiểm là **nguyên nhân** của sự giống nhau. Hai cách giải thích cạnh tranh:

- **H1/H2:** cùng tổ tiên. Khung chung là thừa kế; gene xây khung được giữ nguyên vì đột biến ở đó thường gây chết.
- **H3 (hướng của chủ dự án):** cùng bản thiết kế được tái dùng có chủ đích, không nhất thiết cùng tổ tiên — giống kỹ sư dùng lại một linh kiện tốt cho nhiều sản phẩm không "họ hàng".

Câu hỏi kiểm được, đúng luật của `CLAUDE.md`: **dữ liệu nào phân biệt được hai bên?** R15 trả lời câu này bằng bốn phép phân biệt (mục 2). Trước đó, mục 1 nói phần nào của trực giác đã được đo trên chính chào mào.

---

## 1. Phần "cùng DNA" đã đo được trên chính chào mào

Số liệu từ `docs/05` và README (kết quả A2, 2026-09-08/09; tái lập bằng `pipeline/p02_core_map/`):

| Đại lượng | Giá trị |
|---|---|
| Phần tử bảo tồn qua 363 loài chim ánh xạ được sang chào mào (ngưỡng phủ 0,8) | 248.461 / 2.947.588 (8,43%) |
| Lõi sau gộp trùng lặp | 61,50 Mb = **6,00%** bộ gene chào mào |
| Sau khi căn thêm đoạn ngắn bằng LAST (B1) | 343.781 vùng, 66,63 Mb = 6,50% |
| Phần vùng mã hóa nằm trong lõi (giao base-level) | **59,8%**, làm giàu 9,96 lần so với rải ngẫu nhiên |
| Cùng quy trình trên bộ *tiến hóa nhanh* | 1.133 / 4.054.556 (0,03%) — chênh 281 lần |

Gene đậm lõi nhất theo mật độ là **cụm Hoxa** (Hoxa3/5/11), đúng nhóm gene quy định trục cơ thể. Đây chính là "loại DNA" chủ dự án đang nghĩ tới, và nó có thật trong dữ liệu của dự án.

Ba điều phải đọc kèm, kế thừa `docs/05` phần D:

1. Bảng mật độ gene **chưa qua bộ lọc** `CDS/locus > 0,9` (mô hình gene cụt thổi phồng mật độ). Con số Hoxa chưa được dùng cho phát biểu công khai cho tới khi B2 và bộ lọc D1 chạy xong.
2. Đối chứng "tiến hóa nhanh" chỉ loại được giả thuyết "công cụ căn khớp bừa"; nó **không** chứng minh vùng bảo tồn có chức năng. Đối chứng âm đúng nghĩa (hoán vị vị trí, B2) chưa làm.
3. "Chung" không phải "giống hệt": 6% là phần không đổi; 94% còn lại là phần hai loài chim có thể khác nhau. Trực giác "cùng DNA" đúng cho lớp LÕI của mô hình ba lớp trong README, không đúng cho toàn bộ bộ gene.

Mục 1 trả lời "có một lõi chung không" (có). Nó **chưa** trả lời "lõi chung đó đến từ tổ tiên chung hay từ bản thiết kế chung". Đó là việc của mục 2.

---

## 2. Bốn phép phân biệt, và dữ liệu nghiêng về đâu

Bảng này là bảng Q7 của R15, mở rộng từng dòng với nhãn.

| Dự đoán | Nếu cùng tổ tiên (H1/H2) | Nếu tái dùng thiết kế | Dữ liệu | Nhãn |
|---|---|---|---|---|
| Cây từ dữ liệu độc lập | Hội tụ về một cây, cao hơn ngẫu nhiên | Không có lý do hội tụ; tái dùng cắt ngang cây | Hội tụ | C |
| Gene hỏng, ERV vô dụng | Chung kiểu lỗi ở họ hàng gần, khác kiểu ở dòng xa | Không có động cơ tái dùng đúng kiểu lỗi | Đúng H1 | C |
| Phép kiểm thống kê hình thức | Một tổ tiên chung thắng trên dữ liệu thật | Chưa có mô hình cạnh tranh nào thắng | Thắng áp đảo | C (số odds: U) |
| Công tắc điều hòa chung, cơ quan hội tụ | Công tắc cổ dùng chung, phần còn lại khác nhau | Tái dùng trọn gói thì cả hai nên giống | Đúng H1 | C |

### 2.1 Gene hỏng dùng chung kiểu lỗi

Đây là phép phân biệt mạnh nhất vì nó hỏi về thứ **vô dụng**: một bản thiết kế tái dùng module tốt, không tái dùng lỗi.

- **GULO** (bước cuối tổng hợp vitamin C): ở linh trưởng haplorhine, giả gene mang chỉ thị mất đoạn bảo tồn xuyên suốt phân bộ, khớp một lần bất hoạt ở tổ tiên chung (Mansueto & Good 2024 — #14, C). Linh trưởng mất ~61 triệu năm trước, chuột lang mất ~14 triệu năm trước, **hai thời điểm, hai cơ chế**, khớp hai sự kiện độc lập (Lachapelle & Drouin 2010 — #15, C; tổng quan Drouin 2011 — #16, C).
- **Ở chim**, ca tương đương: khả năng tổng hợp vitamin C ở chim tương quan với phát sinh loài; trong bộ Sẻ, một số loài tổng hợp ở thận và gan, một số chỉ ở gan, một số **mất hẳn** (Chaudhuri & Chatterjee 1969 — #17, C; nhắc lại độc lập ở Liu 2013 — #18, C). Bằng chứng này **yếu hơn linh trưởng**: chỉ ở mức hoạt tính enzyme, chưa có giải trình tự định vị đột biến ở một loài Passeriformes cụ thể (#19, U — khoảng trống thật, đã tìm kỹ).
- **Gene răng chim** (kế thừa R12): cả 6 gene khoáng hóa răng (ENAM, AMBN, AMELX, AMTN, DSPP, MMP20) mang đột biến bất hoạt dùng chung, còn đọc được như giả gene, ở hầu hết các bộ chim còn sống — khớp **một** lần mất răng ở tổ tiên chim hiện đại ~116 triệu năm trước (Meredith 2014 — #21, C, nhãn giữ của R12). Đây là ca "lỗi chung" trên chính lớp chim, và là ca chưa ai kiểm trên bộ Sẻ (xem mục 6).

### 2.2 Thứ bậc lồng nhau từ dữ liệu độc lập

- So khớp 181 cây phân tử với 49 cây hình thái: bất khớp lớn hơn 0 thật, nhưng thấp hơn nhiều so với ngẫu nhiên; cây hình thái và cây phân tử đều xấp xỉ **một** cây nền (Pisani 2007; Patterson 1993 — #7, C).
- Sắc thái bắt buộc: meta-phân tích 2023 trên 32 bộ dữ liệu cho thấy bất khớp topo là "phổ biến" — khớp hoàn hảo không xảy ra (Keating 2023 — #8, C). Đây là đối trọng cho #7, không mâu thuẫn kết luận.
- **SINE** là bằng chứng gần như không homoplasy vì cơ chế chèn độc lập và không đảo ngược (Shedlock & Okada 2000 — #9, C). Ca kinh điển: cá voi lồng vào trong Artiodactyla, gần hà mã nhất, hai nhóm/hai phương pháp cùng kết luận (Shimamura 1999; Nomura & Yasue 1999 — #10, C).
- ERV cũng là công cụ phát sinh loài chuẩn (Johnson & Coffin 1999 — #23, C). Nhưng **không tồn tại** công thức hay con số xác suất đã công bố cho "hai dòng độc lập cùng nhận một chèn ERV tại đúng một tọa độ" — chỉ có lập luận định tính (#24, U — khoảng trống thật, đã tìm bằng PubMed, Consensus, WebSearch, Perplexity).

### 2.3 Phép kiểm thống kê hình thức

- Theobald 2010 (*Nature*) là phép kiểm chọn mô hình duy nhất từng công bố so trực tiếp "một tổ tiên chung phổ quát" với dải mô hình "nhiều tổ tiên độc lập" trên trình tự protein thật, không giả định trước tương đồng trình tự = họ hàng; kết quả nghiêng áp đảo về một tổ tiên chung (#11, C).
- Con số odds cụ thể **không tự xác minh được** qua toàn văn; nguồn thứ cấp mâu thuẫn nhau (#12, U). Không dùng con số nào.
- Có phản biện xuất bản thật (Yonezawa & Hasegawa 2010) và bài trả lời của Theobald 2011 cho một phê bình khác (#13, C). Đây là tranh luận khoa học đang mở về **phương pháp**, không phải về kết luận định tính.

### 2.4 Công tắc điều hòa chung, cơ quan hội tụ

- *eyeless*/Pax6: biểu hiện lạc chỗ tạo mắt hình thái bình thường ở ruồi giấm (Halder 1995 — #1, S cho đúng thí nghiệm gốc); Pax-6 của mực ống cũng tạo mắt lạc chỗ ở ruồi (#2, S); *eyeless* ruồi tiêm vào phôi ếch cũng tạo mắt (#3, S). Gộp ba thí nghiệm, kết luận "Pax6 là công tắc cổ dùng chung xuyên ngành" đạt CONFIRMED với giới hạn ghi rõ: ba nhóm tác giả chồng lấn một phần (#4).
- Nhưng mắt mực và mắt động vật có xương sống **không đồng nguồn**: võng mạc, thủy tinh thể hình thành từ hai lộ trình phôi khác hẳn (Harris 1997 — #5, C); thủy tinh thể chân đầu dùng họ crystallin riêng (Ryu 2023; Sugimoto 2026 — #6, C).
- Đối chiếu (#40, S — lập luận tổng hợp của R15): nếu bộ công cụ được **tái dùng có chủ đích** để tạo mắt ở các dòng không họ hàng, không có lý do để phần còn lại của cơ quan phải xây theo hai lộ trình hoàn toàn khác. Trong khi đó, "công tắc cổ có trước cơ quan + cơ quan hội tụ độc lập" là đúng dự đoán của tiến hóa.

---

## 3. Lập luận phía "thiết kế", trình bày đúng

R15 được giao trình bày trung thực, không dựng bù nhìn. Tóm lại theo #35–#39:

- **Lập luận mạnh nhất:** hội tụ chức năng (mắt camera ở mực và ở người) chứng minh cấu trúc tương tự **có thể** phát sinh không cần tổ tiên chung gần. Vậy tương đồng tự nó không bắt buộc suy ra tổ tiên chung; có thể là một bộ công cụ tốt được tái dùng.
- **Nỗ lực định lượng nghiêm túc nhất:** baraminology thống kê (BDC) — phân cụm thật trên dữ liệu thật, cho phép lồng nhau *trong* một "kind" nhưng không bắt buộc lồng nhau *xuyên* "kind" (#35b).
- **Kênh xuất bản:** các bài này ở *Answers Research Journal*, *Proceedings of the International Conference on Creationism* — không phải bình duyệt khoa học chuẩn; theo luật của `_COMMON.md`, chỉ dùng để mô tả lập luận, không dùng làm chỗ dựa cho claim sinh học (#36).
- **Tranh luận nội bộ thật:** chính tác giả BDC thừa nhận các phương pháp phân cụm khác nhau cho kết quả không nhất quán ở "một số lượng lớn" nghiên cứu (Wood 2021 — #37, C).
- **Phản biện dòng chính** (Sober 2007 — #38, C): cả hai phê bình phổ biến với thiết kế thông minh (không thể kiểm chứng; bị bác bởi thích nghi không hoàn hảo) đều có điểm yếu riêng; dù vậy, khi phân tích cẩn trọng, thiết kế thông minh "không tạo thành lựa chọn thay thế nghiêm túc" cho tiến hóa.
- Phía thiết kế không đồng nhất; young-earth creationism và intelligent design là các lập trường khác nhau (#39, S) — không gộp.

Điểm quyết định (#41, S): mô hình "tái dùng có chủ đích" không có động cơ giải thích vì sao **gene hỏng** lại được tái dùng với đúng kiểu lỗi ở họ hàng gần nhưng khác kiểu ở dòng xa, đúng theo hình cây đã biết từ dữ liệu khác.

---

## 4. "Thiết kế kém": một ca sạch, một ca phức tạp hơn tưởng

- **Dây thần kinh quặt ngược thanh quản:** đo trực tiếp trên hươu cao cổ (Harrison 1980, 1981 — #27, C). Cơ chế phát triển: dây thần kinh hình thành trước khi cổ kéo dài, đường vòng là đặc điểm chung của mọi tetrapoda còn sống; ở khủng long cổ dài, neuron ước tính ≥28 m (Wedel 2012 — #28, C). Chính bài bình duyệt dùng khung "thiết kế không thông minh" (#29, C). Ca này **sạch**.
- **Võng mạc lộn ngược:** giải phẫu "lộn ngược" là thật (Franze 2007 — #30, C). Nhưng khoa học dòng chính tìm ra **hai** lợi ích chức năng độc lập: tế bào Müller dẫn sáng như sợi quang (Franze 2007; Reichenbach & Bringmann 2013 — #31, C) và tiết kiệm không gian ở mắt nhỏ (Kröger & Biehlmaier 2009 — #32, C, lab và cơ chế khác hẳn #31). Luận điểm "kém" **yếu đi nhưng không bị xóa**: điểm mù không được giải thích bởi cả hai, và không có ở mắt chân đầu (#33, C).

Bài học cho dự án: không dùng ví dụ "thiết kế kém" làm luận cứ chính; dùng mục 2.1 và 2.3, là những phép kiểm không phụ thuộc vào đánh giá "tốt/kém".

---

## 5. Nối với R10, R12, R13 và dữ liệu của dự án

Ba brief độc lập, ba góc tiếp cận, cùng một kết quả:

| Brief | Hỏi gì | Trả lời |
|---|---|---|
| R10 / tổng hợp 07 | Có DNA ngoại lai trong chim không? | Có, qua retrovirus, nguồn gốc truy được; quét thành phần của dự án: 0/14.003 vùng mã hóa lệch cả ba thước đo |
| R13 | Bộ gene được thiết kế trông ra sao, đo bằng gì? | Có đối chứng dương thật (watermark JCVI-syn1.0, tái mã hóa Sc2.0) và khung phép đo a–g; **chưa phép đo nào chạy trên dữ liệu chim** |
| R15 (bản này) | Khung chung do tổ tiên hay do thiết kế? | Bốn phép phân biệt đều nghiêng về tổ tiên chung; phía thiết kế chưa nêu dự đoán riêng |

Phép đo gần nhất với khung R13 mà dự án đã chạy là quét thành phần A3 (GC/GC3, cách dùng codon, phổ 4-mer) — âm tính. Phép đo (c) của R13 (pseudogene mục nát khớp họ hàng gần) chính là B5 ở mục 6, chưa chạy.

---

## 6. Verdict cho sổ giả thuyết, và điều kiện để ô H3 được điền

Trong phạm vi R15: dữ liệu nghiêng **rất mạnh** về H1/H2. Phía thiết kế có lập luận thật (mục 3) nhưng **chưa đưa ra dự đoán riêng nào** — khác với H1/H2, kiểm được bằng code — từng thắng trên dữ liệu thật. Ô "dự đoán riêng" của H3 trong `docs/00` Q4 vẫn trống; R15 không điền hộ, chỉ xác nhận nó vẫn trống sau khi đã tìm.

Không viết "đã chứng minh". Viết đúng: *mọi quan sát chủ dự án nêu (khung chung, gene toolkit bảo tồn, ổn định hình thái) đều được H1/H2 giải thích; chưa quan sát nào cần tới H3.*

Điều kiện để một dự đoán H3 được đưa vào bảng: (1) H1/H2 không có dự đoán đó; (2) tính được từ dữ liệu dự án; (3) thắng một phép kiểm chọn mô hình trên dữ liệu thật. Pipeline sẵn sàng chạy khi có.

---

## 7. Việc kiểm được trên chính chào mào

Thứ tự theo chi phí và theo mức trực tiếp với câu hỏi:

1. **B5 — gene răng trên bộ gene chào mào** (`docs/05` B5; phương pháp chuẩn ở tổng hợp 07 mục 5: TBLASTN với ortholog thú/cá sấu làm truy vấn, kiểm thủ công dịch khung và codon dừng). Nếu tìm thấy mảnh vỡ mang cùng kiểu hỏng như các bộ chim khác, đó là phép kiểm 2.1 chạy trên đúng loài của dự án. Chưa ai làm ở bộ Sẻ.
2. **GULOP ở chào mào** — cùng công cụ với B5, lấp đúng khoảng trống #19 của R15 (chào mào thuộc Passeriformes). Dự đoán H1: đột biến bất hoạt giống hoặc lồng được với loài sẻ họ hàng gần; khác hẳn với kiểu hỏng ở linh trưởng và chuột lang.
3. **Lọc D1 và B2** trước khi dùng con số Hoxa ở mục 1 cho bất kỳ phát biểu nào.
4. Vòng 2 của R15 (mục 5 của file thô): số odds Theobald từ bản tác giả tự lưu; công thức xác suất chèn ERV trong giáo trình phát sinh loài; mở toàn văn Harris 1997, Pisani 2007, Patterson 1993, Keating 2023, Sober 2007, Ross 2005 (DOI đang ghi "chưa kiểm").

---

## 8. Giới hạn của bản tổng hợp này

- Đây là tổng hợp **văn liệu**, không phải số đo mới. Số đo chim của dự án (mục 1) vẫn chưa có kiểm định nền và chưa qua bộ lọc D1.
- Bản này **không mở lại nguồn nào**. Nhãn thuần kế thừa từ R15. R15 do subagent Claude viết và **chưa qua review độc lập** (chưa có V-file cho R15, khác A1 đã có V2). Theo bài học 2026-09-08 trong `tasks/lessons.md`, trước khi bất kỳ câu nào ở đây vào bài công bố, cần một lượt mở từng URL.
- Nhiều nguồn trong R15 ghi "DOI chưa kiểm" (qua Consensus): Harris 1997, Ryu 2023, Sugimoto 2026, Pisani 2007, Patterson 1993, Keating 2023, Johnson & Coffin 1999, Sober 2007, Ross 2005, các bài baraminology. Không dùng chúng làm nguồn duy nhất cho claim công khai.
- `research/validate_raw.py` báo NEEDS_FIX cho R15 vì 4 dòng ghi chú (#20, #22, #35, #36) không mang nhãn. Đó là ghi chú kỹ thuật, không phải claim sinh học, nhưng validator không phân biệt. Không sửa tay file thô (luật `CLAUDE.md`); cần quyết định của chủ dự án: thêm nhãn kiểu "GHI CHÚ" vào validator, hoặc giao subagent sửa R15.
- R15 §4 dòng 22 ghi brief nhắc "R14" cho gene răng nhưng nội dung nằm ở R12; R14 (hằng số đột biến) được viết sau. Bản này dùng R12.
