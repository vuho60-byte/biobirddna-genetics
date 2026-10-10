# 12 — Kết quả chạy thật B11 (vị trí câm Hox) và B5 (gene răng) trên bộ gene chào mào

Ngày: 2026-10-10. Người viết: Main. Dữ liệu: `data/b11/` (B11), `data/b5/` (B5 vòng 1 + tóm tắt vòng 3–4),
chạy trên PC chủ dự án (Docker `quay.io/biocontainers/blast:2.16.0--h66d330f_5`, Python 3.14, worktree
`.wt-claude` của nhánh `claude/exciting-ptolemy-0v60q5`). Mã: `pipeline/p08_silent_sites/`, `pipeline/p07_tooth_genes/`
(brief X10, X9). Mốc so sánh văn liệu lấy từ `research/raw/R17` (không mở lại nguồn nào trong bản này).
**Nhãn:** B5 và đếm Hox — STRONG_INFERENCE (một bộ gene mức scaffold, CDS dự đoán). B11 — sau khi lặp trên gà (GRCg7b)
và sẻ vằn (bTaeGut1.4.pri) cùng kết quả: **CONFIRMED trên dữ liệu dự án** (ba bộ gene độc lập, một phương pháp). Không
viết "đã chứng minh".

## 0. Câu hỏi và dự đoán viết trước khi nhìn số

Chủ dự án (2026-10-10): người và chim dùng chung cụm Hox và bộ gene công cụ phát triển → có mở đường cho
"cấu tạo lai tạo từ DNA xuyên loài" không? Dự án không trả lời bằng ý kiến mà bằng hai phép đo có dự đoán viết trước:

| Phép đo | Dự đoán nếu thừa kế từ tổ tiên chung (H1/H2) | Dự đoán nếu chép/ghép gần đây |
|---|---|---|
| B11: độ giống vị trí câm (base 3 codon 4-fold) người–chào mào ở Hox, so với gene đối chứng | tiến về mức ngẫu nhiên tính từ thành phần base; Hox giống đối chứng | gần 1,0; Hox lệch khỏi đối chứng |
| B5: mảnh gene răng hỏng ở chào mào mang tổn thương khớp các bộ chim khác | có, nếu công cụ đủ nhạy (đối chứng dương phải đạt) | không có lý do mang gene hỏng cùng kiểu lỗi |

## 1. B11 — vị trí câm: Hox không khác gì gene đối chứng, không gene nào "giống hệt"

Đầu vào: 60/73 gene người có CDS RefSeq (39 HOX có thật của người + 11 gene công cụ + 10 gene đối chứng; 13 tên
HOX không tồn tại ở người đúng như kỳ vọng). 41 gene có CDS chào mào cùng tên trong chú giải D01; 39 cặp có ≥ 20
vị trí 4-fold trong vùng lõi (cửa sổ 21 aa ≥ 50 % giống).

| Lớp | Số gene đo được | aa giống (lõi), trung vị | Giống vị trí 4-fold, trung vị | Mức ngẫu nhiên, trung vị | Nhãn |
|---|---|---|---|---|---|
| Hox | 22 | 85,3 % | 48,9 % | 30,7 % | STRONG_INFERENCE |
| Gene công cụ (PAX6, SOX9, RUNX2, SHH, BMP4, TBX4/5, PITX1, FGF8, GLI3, HAND2) | 11 | 88,7 % | 53,4 % | 32,5 % | STRONG_INFERENCE |
| Đối chứng housekeeping (ACTB, GAPDH, ENO1, LDHA, RPLP0, TPI1) | 6 | 91,7 % | 49,2 % | 27,5 % | STRONG_INFERENCE |

Đọc:
- **Không gene nào trong 39 gene có độ giống vị trí câm ≥ 0,9.** Cao nhất là HOXD8 74,4 % (43 vị trí) và HOXA11 65,1 %;
  gene đối chứng ACTB, protein giống 100 %, chỉ 59,0 % ở vị trí câm. Kịch bản "chép/ghép gần đây" (vị trí câm ≈ 1,0)
  **không khớp** bất kỳ gene nào, kể cả Hox.
- **Hox = đối chứng**: trung vị 48,9 % so với 49,2 %. Nếu Hox được ghép từ nguồn khác, Hox phải lệch khỏi phần còn
  lại của bộ gene ở vị trí câm; nó không lệch.
- Mọi lớp đều **trên mức ngẫu nhiên ~30 %** nhưng xa 1,0: đúng dạng "thừa kế lâu đời + ràng buộc codon", khoảng
  cách JC 0,5–1,5 thay thế/vị trí (nhiều gene gần bão hòa). R17 #40–#43 cho mốc văn liệu: Ks người–chuột 0,6–0,73,
  Ks người–gà trung vị 1,66 bão hòa (STRONG_INFERENCE); người–chào mào kỳ vọng ≈ người–gà — số đo ở đây cùng bậc,
  thấp hơn một phần vì chỉ tính cột lõi và JC đánh giá thấp khi bão hòa.
- Phát hiện phụ: 19/39 HOX người **không có CDS chào mào cùng tên** (HOXA2, A10, A13, B1, B5, B9, C4–C6, C9–C11, C13,
  D1, D3…) và 7 gene có CDS chào mào **cụt** (156–333 nt: HOXA4, B4, B7, B13, C8, C12, D13). Đây là lỗ hổng chú giải
  của bộ gene mức scaffold, khớp cảnh báo D1 (mô hình gene cụt) ở `docs/05` — không phải mất gene thật. Đếm số Hox thật
  bằng tblastn 39 protein HOX người lên bộ gene (vòng 7, lệnh con `census`: mỗi vùng genome gán cho truy vấn tổng
  bitscore cao nhất; `data/b5/r7_hox_census_summary.tsv`): **31/39 HOX người có vùng mạnh không mơ hồ** trên bộ gene
  chào mào; HOXA6 nằm chung vùng với HOXA5 (hai gene cách < 3 kb, vùng 339–345 kb của VWYP01010181.1 gán cho HOXA5 với
  HOXA6 đứng nhì); HOXB7, HOXC8 có vùng mạnh nhưng cụt ở mép scaffold; chỉ HOXB1, HOXC4, HOXC5, HOXD1 còn yếu (bitscore
  60–81). Bốn cụm nằm đúng chỗ: HOXA trên VWYP01010181.1 (284–392 kb, 10 gene thứ tự A13→A1), HOXD trên VWYP01016378.1
  (7–64 kb, 9 gene D4→D13), HOXB và HOXC rải trên nhiều scaffold ngắn. Kết luận: chào mào có ít nhất 34/39 gene Hox của
  người, thứ tự trong cụm giữ nguyên; 19 ô "không có CDS" của B11 là lỗ hổng chú giải, không phải mất gene.
  STRONG_INFERENCE (một bộ gene; ô trống của người như HOXC5 không kiểm được bằng truy vấn người).

**Lặp trên gà và sẻ vằn (cùng CDS người, cùng tham số; `data/b11/chicken/`, `data/b11/zebrafinch/`):**

| Bộ gene | Gene đo được | Hox: giống 4-fold / ngẫu nhiên | Công cụ: giống 4-fold / ngẫu nhiên | Đối chứng: giống 4-fold / ngẫu nhiên | Gene ≥ 0,9 |
|---|---|---|---|---|---|
| Chào mào GCA_013400435.1 (scaffold) | 39 | 48,9 % / 30,7 % (22 gene) | 53,4 % / 32,5 % (11) | 49,2 % / 27,5 % (6) | 0 |
| Gà GRCg7b (nhiễm sắc thể) | 54 | 49,3 % / 32,6 % (34) | 52,4 % / 31,8 % (11) | 52,6 % / 28,9 % (9) | 0 |
| Sẻ vằn bTaeGut1.4.pri (nhiễm sắc thể) | 55 | 49,8 % / 33,5 % (36) | 52,8 % / 32,6 % (11) | 51,0 % / 28,3 % (8) | 0 |

Ba bộ gene, ba phòng thí nghiệm lắp ráp khác nhau, cùng bức tranh: Hox ở vị trí câm giống người 49–50 %, bằng gene
đối chứng, cách xa 1,0, trên mức ngẫu nhiên ~30 %. Bộ gene gà và sẻ vằn chú giải đầy đủ hơn nên đo được 34–36 Hox
thay vì 22.

Verdict B11: dữ liệu **loại kịch bản chép/ghép gần đây** cho toàn bộ gene đo được ở cả ba loài chim và cho thấy Hox
không khác phần còn lại của bộ gene ở vị trí câm. **CONFIRMED trên dữ liệu dự án** (ba bộ gene độc lập; vẫn là một
phương pháp — phép kiểm độc lập thứ hai là cây gene so cây loài, chưa làm).

## 2. B5 — gene răng: đối chứng MEPE không đạt, năm gene men/ngà răng không hit, MMP20 chỉ hit chéo họ MMP

Truy vấn: protein cá sấu Mỹ RefSeq ENAM (1092 aa), AMBN (408), AMELX (199), AMTN (328), DSPP (797), MMP20 (485);
đối chứng MEPE (847), SPARC (298); GULO **không có bản ghi tên GULO** ở cá sấu lẫn gà trong RefSeq (vòng 4) nên vòng 5
dùng Gulo chuột cống.

| Mục | Kết quả | Đọc |
|---|---|---|
| SPARC (đối chứng) | 1 locus, phủ 0,95, 95 % id, 0 codon dừng, 0 lệch khung; GFF: đúng gene *Sparc* | ĐẠT — quy trình thấy được gene nguyên vẹn |
| MEPE (đối chứng, họ SCPP như gene răng) | **không tìm thấy** ở e-value 1e-5 | KHÔNG ĐẠT — độ nhạy cho họ SCPP (protein giàu P/Q/S, tiến hóa nhanh) không đủ với truy vấn cá sấu ~250 Ma |
| ENAM, AMBN, AMELX, AMTN, DSPP | không hit | **không kết luận** — khi đối chứng cùng họ không thấy, âm tính không nói gì về giả gene |
| MMP20 | 11 locus, 38–70 % id; GFF: **mọi locus chồng gene Mmp họ hàng đã chú giải** (Mmp2, 7, 8, 9, 10, 11, 15, 16, 17, 24, 27, 28, 3) trừ một **giả gene chưa tên** VWYP01031341.1:548624–554768 trong cụm MMP | hit chéo paralog; 2 "codon dừng có chỗ dựa" của báo cáo nằm trong vùng *Mmp27* và trong giả gene đó — **không phải tổn thương MMP20** |
| Giả gene chưa tên (vòng 3, bộ MMP người) | MMP13 8 HSP, tổng 740 bit, 89 % id, phủ aa 41–471; MMP20 1 HSP, 55 bit, 51 % | giả gene là **tàn dư MMP13**, không phải MMP20 (`data/b5/r3_mmp_family_summary.txt`) |
| Tìm theo vị trí láng giềng (vòng 6, e-value 1e-2, mốc SPARCL1/SPP1/IBSP/DMP1) | mốc tìm thấy đúng cụm: DMP1 4,2 kb, IBSP 15,3 kb và một hit MEPE yếu (25 % id) 23,7 kb trên cùng scaffold VWYP01026701.1 dài 38,7 kb; SPARCL1 ở mép scaffold khác 36 kb; ENAM/AMBN/AMELX/AMTN/DSPP 0 HSP | vùng DSPP rơi vào khoảng trống giữa hai scaffold ngắn; cụm men răng không có hit → **không kết luận**, do lắp ráp (`data/b5/r6_scpp_synteny_summary.txt`) |
| GULO (vòng 5, Gulo chuột cống NP_071556.2 vì cá sấu/gà không có bản ghi tên GULO) | 6 mảnh trên 6 scaffold khác nhau, mỗi mảnh 1 HSP phủ 7–14 % protein, 56–79 % id, tổng phủ ~0,62; nhiều mảnh sát đầu scaffold; 0 codon dừng có chỗ dựa | vùng GULO **vỡ trong lắp ráp**; còn chức năng hay giả gene: **không kết luận** được từ bộ gene này (`data/b5/r5_gulo_summary.txt`); khoảng trống #19 của R15 vẫn mở |

Đọc theo thứ tự khả năng của `report.md`: (1) không có tổn thương nào của gene mục tiêu để xét; (2) hai `*` "có chỗ
dựa" là giả tạo do HSP paralog mở rộng — lọc cánh 10 aa không đủ khi paralog giống 40–60 %; (3) ranh exon: 11 cặp
`intron_like` đều là exon của gene Mmp chức năng; (4) **không có mảnh MMP20 phát hiện được** trên bộ gene chào mào ở
ngưỡng này. "Không tìm thấy ≠ không tồn tại": Meredith 2014 (R12 #23, CONFIRMED) báo tàn dư MMP20 đọc được ở hầu hết
bộ chim; ở chào mào nó có thể nằm trong khoảng trống lắp ráp (N50 ~218 kb; vùng MMP7–MMP27 chỉ cách 43 kb không có chú
giải nào) hoặc mục nát dưới ngưỡng tblastn với truy vấn cá sấu. Phép thử B5 **chưa phân biệt được gì** cho H1/H2/H3.

Bài học công cụ (vào vòng sau của p07): (a) bước GFF-overlap + kiểm paralog (best reciprocal hit) phải là bước
chuẩn trước khi gọi "tổn thương"; (b) truy vấn nên là **tàn dư của chính chim** (tọa độ từ phụ lục Meredith 2014 ở gà/sẻ
vằn, hoặc tblastx bộ gene–bộ gene) thay cho protein cá sấu; (c) đối chứng cùng họ (MEPE) đã làm đúng việc của nó: chặn
diễn giải quá đà.

## 3. R17 — "cùng sai số" áp vào người: sáu dự đoán, cả sáu nghiêng về dòng linh trưởng

`research/raw/R17-nguoi-cung-sai-so-voi-linh-truong.md` (49 dòng; 27 C / 23 S / 12 U; 34 bài bình duyệt đã mở; 8 trang
phía thiết kế trình bày nguyên văn). Tóm tắt mục 4 câu 7: lỗi gene vô dụng chung chỗ trong nhóm, khác chỗ ngoài nhóm
(GULOP 59 thay thế nút Haplorhini; urate oxidase CGA→TGA chung người + vượn lớn, gibbon hỏng khác — CONFIRMED); chèn
trung tính trùng tọa độ hàng nghìn, nhiễm độc lập chỉ trùng 0,4–4,2 % (CONFIRMED / STRONG_INFERENCE); lỗi riêng người
ít và trẻ (CMAH, MYH16 — CONFIRMED); vị trí câm tăng theo thời gian, dS người–tinh tinh 1,2–1,6 %, người–gà bão hòa
(CONFIRMED / STRONG_INFERENCE); NST 2 đủ ba dấu vết (CONFIRMED). Ô "dự đoán riêng của phía 'dòng tạo riêng' đã kiểm và
thắng" **vẫn trống** (UNVERIFIED — ô trống thật). R17 cũng sửa một kỳ vọng từ trí nhớ của Main: Ks người–chuột là
0,6–0,73, không phải 0,3–0,5.

## 4. Trả lời cho câu hỏi của chủ dự án

Năm dấu vết mà một bộ gene "ghép xuyên loài" phải để lại, trạng thái sau hôm nay:

| Dấu vết | Phép đo | Kết quả | Nhãn |
|---|---|---|---|
| Đoạn ghép mang chữ ký thành phần lạ | A3 quét GC/codon/4-mer trên 14.003 CDS | 0 đoạn lệch chuẩn | CONFIRMED trên dữ liệu dự án |
| Vị trí câm giống hệt nguồn | B11, người so chào mào, gà, sẻ vằn (39–55 gene mỗi loài) | không gene nào ≥ 0,9 ở cả ba loài; Hox = đối chứng | CONFIRMED trên dữ liệu dự án |
| Không sao chép gene hỏng cùng kiểu lỗi | R15 (chim), R17 (người) văn liệu; B5 chào mào | văn liệu: lỗi chung lồng theo cây; B5: chưa đo được (độ nhạy) | CONFIRMED (văn liệu) / chưa kết luận (B5) |
| Vùng khác nhau cho cây khác nhau | so cây từng gene với cây loài | chưa làm | — |
| Cụm ghép không giữ trật tự gene/intron | so GFF chào mào với người quanh bốn cụm Hox | chưa làm (rẻ) | — |

Hai dấu vết đo được đều không ủng hộ "ghép"; hai dấu vết văn liệu nghiêng hẳn về tổ tiên chung; hai dấu vết chưa đo.
Giả thuyết "lai tạo xuyên loài" và "người là dòng tạo riêng" hiện **chưa có dự đoán riêng nào thắng trên dữ liệu thật**
— cùng tình trạng với ô H3 trong `docs/00` Q4. Điều kiện để lật: một gene hoặc một vùng có vị trí câm gần giống hệt
giữa hai loài xa (B11 lặp trên gà/sẻ vằn vẫn không thấy), hoặc cây gene mâu thuẫn cây loài vượt mức ngẫu nhiên có kiểm định.

## 5. Việc tiếp theo

1. ~~B11 lặp trên gà và sẻ vằn~~ — đã làm, cùng bức tranh. Phép kiểm độc lập thứ hai: cây từng gene so cây loài.
2. **B5 vòng 2**: truy vấn bằng tàn dư gene răng của chim (phụ lục Meredith 2014) hoặc tblastx; thêm bước GFF-overlap +
   kiểm paralog vào p07; thử e-value nới (1e-2) có kiểm đối chứng MEPE; GULO cần bộ gene sẻ lắp ráp tốt hơn (6 exon tìm được nằm trên 6 scaffold rời).
3. ~~Đếm Hox thật~~ — đã làm (34/39, xem mục 1). Còn mở: ô trống của lưới Hox người (A8, A12, B10–12, C1–3, C7, D2, D5–7) có trống ở chim không — cần truy vấn từ loài có gene đó (cá).
4. **R17 vòng 2** theo 12 câu ở mục 5 file thô; V-file review độc lập cho R15 và R17.
5. Ghi vào `docs/00` Q4: thêm dòng B11 vào bảng dự đoán phân biệt (sửa vòng sau, hậu tố `*(sửa theo B11)*`).

## 6. Giới hạn của bản này

- Toàn bộ số đo đến từ **một** bộ gene mức scaffold với chú giải tự động; CDS cụt và thiếu làm giảm số gene đo được,
  không làm lệch chiều kết luận (vị trí câm của gene cụt vẫn không gần 1,0).
- B5 là kết quả **âm tính không kết luận**; không được dùng làm bằng chứng cho hay chống bất kỳ H nào.
- `report.md` trong `data/b11` và `data/b5` được tạo lại trong kho từ TSV đã kiểm sha256 khớp PC (chuỗi đường dẫn khác
  nên sha báo cáo khác bản PC); `human_cds.fna`, `alignments.txt`, `queries.faa`, `tblastn.tsv` giữ trên PC, sha256 ghi
  trong `data/registry.csv` (D10-b11, D09-b5) và ghi chú registry.
- Chưa ai review độc lập mã p07/p08 ngoài unittest; vòng 3–5 của B5 chỉ có tóm tắt lọc (`data/b5/r3_*.txt`, `r5_*.txt`), dữ liệu đầy đủ ở `data/r3–r5` trên PC.
