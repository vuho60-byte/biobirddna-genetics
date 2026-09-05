# Ranh giới khoa học, an toàn và ngôn ngữ công bố

Đọc file này trước mọi thay đổi model hoặc report.

---

## 1. Bốn cấp không được trộn

| Cấp | Dữ liệu | Được phép nói |
| --- | --- | --- |
| Quan sát | Ảnh/video/feature | “Có vùng trắng/nâu/xám/đỏ thay đổi…” |
| Kiểu hình | T02 + timeline | “Kiểu hình gợi ý MEL-PGR/MEL-LEU…” |
| Mô hình di truyền | Phả hệ + toàn bộ đời con | “Phù hợp hơn với lặn thường/liên kết Z…” |
| Cơ chế phân tử | DNA + đồng phân ly + tái lập | “Locus/biến thể được hỗ trợ…” |

Không nhảy cấp.

---

## 2. Câu bị cấm từ ảnh

- “Con này mang gene lặn.”
- “Đây là gene TYR/SLC45A2/PMEL.”
- “Màu này liên kết giới tính.”
- “Ghép với con X sẽ ra 50%.”
- “Albino đã xác nhận.”
- “Thiếu chất gây màu này.”

Ảnh có thể đo phenotype, không đọc genotype.

---

## 3. Câu được phép

- “Ảnh Q3 cho thấy khoảng 18% vùng lông sáng/trắng.”
- “Diễn biến tăng qua hai mùa lông làm MEL-PGR được ưu tiên hơn MEL-LEU.”
- “Trong các mô hình đã kiểm tra, dữ liệu đời con hiện phù hợp hơn với lặn liên kết Z.”
- “Kết quả này là mô hình ứng viên G2; cần phép lai độc lập và giới tính đời con để nâng bằng chứng.”
- “Nếu giả định cha N/m và mẹ m/m theo mô hình lặn thường, xác suất kỳ vọng phenotype mục tiêu là 50%.”

---

## 4. Evidence grade và claim ceiling

| Grade | Claim ceiling |
| --- | --- |
| G0 | Không kết luận di truyền |
| G1 | Mô hình thăm dò |
| G2 | Mô hình ứng viên |
| G3 | Mô hình phân ly được hỗ trợ |
| G4 | Mô hình mạnh, chưa xác nhận locus |
| G5 | Locus/biến thể được dữ liệu hỗ trợ |
| G6 | Cơ chế xác nhận mạnh |

Numerical posterior không tự nâng grade.

Ranh giới này được engine cưỡng chế bằng code, không chỉ bằng quy ước: mọi
kết quả đi qua `dispatch()` (Python `service.py`, Node `index.mjs`) đều bị ép
qua evidence boundary (`evidence.py`) — claim mạnh nhất được phép của từng
operation được tính bằng `safe_claim_language` và gắn vào
`limits.claim_language`; claim vượt trần grade làm `assert_claim_allowed` ném
`EvidenceBoundaryError` và không có output nào được phát ra. Trần hiện tại:
`cross`/`rank-mates`/`simulate-recovery` = `simulate_cross` (G0), `pigment` =
`phenotype_hypothesis` (G0), `normalize-records` = `observation_record` (G0),
`infer` = `candidate_inheritance_model` (G2).

Feedback artifact cũng không nâng evidence grade và không biến metric vận hành
thành bằng chứng sinh học. Health gate 3.0 của AutoResearch không phải score
gene/model BioBirdDNA; adapter chỉ được giữ trạng thái chờ reviewer.

---

## 5. Bias phải ghi nhận

- Chỉ lưu con màu lạ.
- Bỏ qua trứng không nở/chết non.
- Chấm phenotype trước thay lông đầu.
- Nhiều ảnh một con thành nhiều mẫu.
- Chọn cặp sau khi đã biết outcome rồi báo như test độc lập.
- Trộn nhiều dòng/cặp.
- Sai giới tính.
- Pedigree do lời kể không kiểm chứng.
- Chỉ theo dõi lứa thành công.
- Seller claim ảnh hưởng reviewer.

---

## 6. Cận huyết và phúc lợi

- Score màu không override block phúc lợi.
- Default chặn trực hệ và full sibling.
- Ngưỡng F còn lại là policy có reviewer, không phải số “an toàn tuyệt đối”.
- Không dùng skill để khuyến khích lặp một trống quá nhiều thế hệ.
- Cần theo dõi fertility, hatchability, survival và dị tật, không chỉ màu.
- Cặp có nguy cơ sức khỏe/phúc lợi phải `welfare_blocked:true` dù xác suất target cao.

---

## 7. Molecular boundary

G5 tối thiểu cần:

- đúng cá thể/mẫu;
- chain of custody;
- genotype QC;
- phenotype QC;
- cosegregation;
- sample size và family structure phù hợp;
- loại trừ lỗi sample swap;
- báo variant/locus cụ thể.

Association một ca không đủ causal.

---

## 8. Pigment chemistry boundary

RGB/CIELAB/spectral appearance không đồng nghĩa trực tiếp với phân tử pigment.

Để nói thành phần hóa học cần phương pháp phù hợp như:

- reflectance spectroscopy;
- Raman;
- HPLC/LC-MS;
- microscopy/chemical assays.

Skill v0.3.0 luôn trả `biochemical_pigment_status:unconfirmed`.

---

## 9. Research pairing vs production pairing

### Research pairing

Mục tiêu: outcome giữa model khác nhau rõ, EIG cao.

### Production pairing

Mục tiêu: phenotype mong muốn, sức khỏe, đa dạng, quản trị dòng.

Không mặc định cặp nghiên cứu tốt là cặp sản xuất tốt.

---

## 10. Publication language template

### G1–G2

> Dữ liệu hiện tại chỉ hỗ trợ một mô hình thăm dò/ứng viên. Tỷ lệ đời con còn nhỏ và chưa đủ để loại trừ các mô hình cạnh tranh.

### G3

> Kết quả từ nhiều phép ghép độc lập phù hợp với mô hình phân ly …, nhưng chưa xác định locus hoặc biến thể phân tử.

### G4

> Phả hệ nhiều thế hệ và test cross hỗ trợ mạnh mô hình …; xác nhận gene vẫn cần dữ liệu phân tử đồng phân ly.

### G5

> Biến thể/locus … đồng phân ly với phenotype trong tập gia đình đã phân tích; cần tái lập độc lập hoặc bằng chứng chức năng để xác nhận cơ chế mạnh hơn.

---

## 11. Safety test bắt buộc

Một patch fail nếu:

- pigment output có `confirmed_gene`;
- inference gọi posterior là gene probability;
- blocked pair vẫn `eligible:true`;
- code chấp nhận Z-linked genotype không có sex;
- report bỏ missing data;
- test gene claim G4 được cho phép;
- full sibling không block theo default.
