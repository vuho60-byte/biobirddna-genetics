# R17 — Người có "cùng sai số" với linh trưởng không? Dữ liệu nào phân biệt "cùng tổ tiên" với "dòng được tạo riêng"?

Đọc `research/briefs/_COMMON.md`. Output: `research/raw/R17-nguoi-cung-sai-so-voi-linh-truong.md`.

**Công cụ trong phiên này (khác _COMMON vì PubMed/Consensus/Perplexity MCP không kết nối được):** `WebSearch`, `WebFetch`, `mcp__firecrawl__firecrawl_search`, `mcp__firecrawl__firecrawl_scrape`, `mcp__tavily_MCP__tavily_search`/`tavily_extract`. Mỗi nguồn dùng cho nhãn CONFIRMED/STRONG_INFERENCE phải được **MỞ** (fetch trang abstract PubMed/PMC/tạp chí), không chỉ dựa vào đoạn tóm tắt của máy tìm kiếm — lesson 2026-09-08 (Antigravity bịa nguồn). Ghi "Perplexity: 0 lượt (không có)" trong search log.

## Câu hỏi thật đằng sau brief này

Chủ dự án vừa rút ra từ R15 một nguyên lý đúng: **chỉ có tương đồng (cùng tổ tiên) mới sinh ra cùng một sai số ở cùng một chỗ**. Anh ấy đặt hai giả thuyết đối lập cho con người: (a) người thuộc dòng linh trưởng như mọi loài khác; (b) người **không** từ dòng linh trưởng mà là một dòng được tạo/ghép riêng ("một thí nghiệm"). Dự án đối xử với cả hai bằng cùng một câu: **dữ liệu nào phân biệt được?** Không kết luận bằng ý kiến. Không chế nhạo.

## Câu hỏi phải trả lời

1. **GULO ở người.** R15 đã có #14–#16 (GULOP giả gene ở haplorhine, indel bảo tồn, 61 Ma vs chuột lang 14 Ma). Bổ sung riêng cho người: đột biến bất hoạt nào **chung giữa người và tinh tinh/khỉ**, có đột biến nào **riêng người** không, và tỉ lệ đột biến chung/riêng khớp cây linh trưởng không? Nêu nguồn giải trình tự (không chỉ tổng quan).
2. **Nhiễm sắc thể 2 của người.** Xác minh: vị trí dung hợp 2q13 có lặp telomere đầu-đối-đầu (TTAGGG)n(CCCTAA)n; tâm động thoái hóa ở 2q21; hai nhiễm sắc thể tương ứng ở tinh tinh/gorilla/orangutan. Nêu công trình gốc (tác giả, năm; tạp chí chỉ ghi khi đã mở nguồn). Nêu cả **phản bác đã công bố** của phía thiết kế (ví dụ lập luận "vùng dung hợp nằm trong gene DDX11L2", "thiếu lặp telomere") và trả lời của giới di truyền học.
3. **Retrovirus nội sinh (ERV) cùng vị trí.** Số lượng ERV ở người có bản sao **cùng tọa độ** ở tinh tinh và các vượn khác là bao nhiêu (theo bộ gene tinh tinh 2005 hoặc nghiên cứu sau)? Có bao nhiêu ERV **riêng người** (HERV-K HML-2) và điều đó nói gì? Xác suất trùng ngẫu nhiên được ai tính? Phía thiết kế lập luận "ERV có chức năng nên không phải tai nạn" — nêu đúng và nêu trả lời (chức năng được chiêu mộ sau không xóa dấu chèn chung vị trí).
4. **Giả gene và yếu tố chèn khác chung vị trí người–vượn.** Một hoặc hai ca có giải trình tự rõ: giả gene đã xử lý (processed pseudogene) hoặc Alu/L1 chèn cùng vị trí ở người và tinh tinh; ca kinh điển nhất được trích dẫn nhiều.
5. **Phép thử vị trí câm.** Nếu bộ gene người được "chép/ghép" gần đây từ một mẫu khác, các vị trí câm (codon position 3, 4-fold degenerate) sẽ gần giống hệt nguồn; nếu thừa kế lâu đời dưới chọn lọc, protein giữ nhưng vị trí câm bão hòa. Tìm con số đã công bố: độ phân kỳ vị trí đồng nghĩa (dS hoặc Ks) người–tinh tinh, người–chuột, người–gà (kỳ vọng thô: người–tinh tinh ~1 %, người–chuột ~0,3–0,5 thay thế/vị trí, người–gà bão hòa — **đây là kỳ vọng từ trí nhớ của Main, subagent phải thay bằng số có nguồn hoặc đánh UNVERIFIED**). Dự án sẽ tự tính người–chào mào ở `pipeline/p08_silent_sites/`; R17 chỉ cần mốc so sánh có nguồn.
6. **Lập luận mạnh nhất của phía "người là dòng tạo riêng".** Trình bày trung thực (tài liệu có bình duyệt nếu có; nếu chỉ có kênh không bình duyệt thì ghi rõ), và nêu dự đoán kiểm được nào họ đưa ra. Không dựng bù nhìn.
7. **Câu quyết định.** Liệt kê các dự đoán **phân biệt được** giữa (a) và (b) ở mục trên; với từng dự đoán ghi dữ liệu hiện có ủng hộ bên nào, nhãn C/S/U. Nếu nghiêng hẳn một bên, nói thẳng, nhưng nêu bên kia lập luận gì.

## Ranh giới bắt buộc
- **Cấm** dùng nguồn không bình duyệt làm chỗ dựa cho claim thực tế (được trích để trình bày lập luận của một phía, ghi rõ "không bình duyệt").
- **Cấm** chế nhạo hay hạ thấp quan điểm nào.
- Nhãn C/S/U cho **mọi** dòng trong bảng phát hiện, kể cả dòng ghi chú kỹ thuật (validator `research/validate_raw.py` sẽ báo NEEDS_FIX nếu một dòng không có nhãn — R15 từng vướng lỗi này).
- Mục 3 (dataset) có thể ghi "không có dataset mới" nhưng vẫn phải có bảng.
- Không sửa file nào ngoài file thô và MỘT dòng `tasks/primordial-dna-handoff.md`.
