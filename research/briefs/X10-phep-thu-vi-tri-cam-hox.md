# X10 — Việc B11: phép thử vị trí câm cho Hox và bộ gene công cụ, người so với chào mào

Chỉ tạo/sửa trong `pipeline/p08_silent_sites/**` (+ 1 dòng `tasks/todo.md`, 1 dòng `tasks/primordial-dna-handoff.md`). **Không có lệnh xóa nào.** Python 3.11 **stdlib**. Không chạy mạng ở bước viết code.

Người thực hiện: Main (phiên 2026-10-10), theo lựa chọn "2 và 3, sau đó 1" của chủ dự án.

## Câu hỏi khoa học

Thảo luận 2026-10-10 (docs/00 Q4; tổng hợp 11): người và chim dùng chung cụm Hox và bộ gene công cụ phát triển (Pax6, Sox9, Runx2, Shh, Bmp4, Tbx5, Pitx1 — đều có mặt và giàu lõi trong `data/a2/gene_enrichment.tsv`). Chủ dự án hỏi: điều đó có mở đường cho giả thuyết "cấu tạo lai tạo từ DNA xuyên loài" không? Dự án không trả lời bằng ý kiến mà bằng một phép đo có **dự đoán viết trước**:

| Kịch bản | Protein | Vị trí câm (base 3 của codon 4-fold) |
|---|---|---|
| Thừa kế từ tổ tiên chung ~320 Ma, chọn lọc giữ protein (H1/H2) | giống cao ở vùng lõi | bão hòa → tiến về mức ngẫu nhiên tính từ thành phần base |
| Chép/ghép gần đây giữa hai bộ gene | bất kỳ | gần giống hệt (~1,0) |

Vị trí 4-fold là nơi **mọi** thay đổi base thứ 3 đều không đổi axit amin, nên chọn lọc lên protein không giữ được nó; chỉ thời gian ngắn (hoặc chép trực tiếp) mới giữ được. Đây là phép thử phân biệt sắc nhất trong năm dấu vết đã liệt kê cho chủ dự án, và tính được ngay từ dữ liệu D01 (`cds_from_genomic.fna.gz`) cộng CDS người tải từ NCBI (vài trăm KB).

## Dữ liệu
- `data/ncbi/GCA_013400435.1/..._cds_from_genomic.fna.gz` — CDS chào mào (D01), header `[gene=Hoxa1] [protein_id=...]`.
- CDS người: RefSeq mRNA qua E-utilities (`esearch` nuccore `GENE[Gene Name] AND "Homo sapiens"[Organism] AND refseq[Filter] AND biomol_mrna[PROP]` → `efetch rettype=fasta_cds_na`), giữ CDS dài nhất đúng `[gene=]`. **Không ghi accession từ trí nhớ**; manifest + 1 dòng `data/registry.csv` (sha256).
- Gene: 39 HOX người (gọi đủ 52 tên A1–D13, tên không có sẽ `not_found`), 11 gene công cụ, 10 gene đối chứng housekeeping (ACTB, GAPDH, EEF1A1, HSPA8, PGK1, TPI1, ENO1, LDHA, ALDOA, RPLP0).

## Việc — `pipeline/p08_silent_sites/silent_sites.py`, ba lệnh con

1. `fetch-human --genes ... --email X --out human_cds.fna --manifest human_manifest.tsv [--registry data/registry.csv]`
2. `compare --human human_cds.fna --bulbul <cds_from_genomic.fna.gz> --out silent_sites.tsv [--alignments alignments.txt]`
   - Khớp gene theo tên (không phân biệt hoa/thường), mỗi bên giữ CDS dài nhất; loại CDS không chia hết 3, < 30 codon, có codon dừng nội (ghi `status`).
   - Dịch protein → căn toàn cục Needleman–Wunsch/Gotoh (BLOSUM62, gap 11/1, tự cài) → ánh xạ về codon.
   - **Cột lõi** = cột không gap nằm trong cửa sổ ±10 aa có ≥ 50 % aa giống (loại vùng căn mơ hồ ngoài homeodomain, vùng lặp poly-Ala/Pro của Hox). Thống kê chính trên cột lõi; thống kê "all" ghi kèm.
   - Mỗi gene: `aa_ident`, `ident_p1/p2/p3`, `n_4fold`, `ident_4fold`, `chance_4fold` (Σ f_người(b)·f_chào mào(b) tại vị trí 4-fold), `jc_4fold` (JC69, ≥ 0,75 → inf), GC3 hai bên; `status=ok` khi ≥ 20 vị trí 4-fold lõi.
3. `report --tsv silent_sites.tsv --out report.md` — bảng theo lớp hox / toolkit / control, trung vị từng lớp, so Hox với đối chứng, verdict từng gene (≥ 0,9 → "gần giống hệt, kiểm lại"; trong ±0,15 quanh ngẫu nhiên → "khớp thừa kế lâu đời"; giữa → "thừa kế + ràng buộc codon"), và **phần cảnh báo bắt buộc** (một cặp loài; CDS dự đoán; GC khác nhau; codon bias; không phân biệt được với kịch bản thiết kế cố tình ngẫu nhiên hóa vị trí câm — kịch bản đó phải tự nêu dự đoán riêng).

## Test (`tests/test_silent_sites.py`, unittest, dữ liệu nhúng)
BLOSUM62 đối xứng, đường chéo dương; căn giống hệt không gap, chèn 3 aa thành đúng 3 gap; bảng mã và 4-fold (Leu CTN đúng, Leu TTR sai, Ile sai, ký tự mơ hồ sai); JC69; CDS giống hệt → mọi độ giống 1,0; đổi base 3 của mọi codon 4-fold → protein giữ nguyên, `ident_4fold` 0, JC inf; ngẫu nhiên hóa base 3 → `ident_4fold` < 0,6, `chance` trong (0,15; 0,6); mặt nạ lõi loại đúng khối 30 aa thay bằng axit amin tương tự; `compare` end-to-end trên FASTA nhúng (isoform dài nhất, codon dừng nội → bỏ, gene thiếu → `no_bulbul_cds`, lớp hox/toolkit/control); `report` in đủ verdict và cảnh báo; `fetch-human` với mạng giả lập (chọn CDS dài nhất đúng gene, bỏ gene sai, `not_found`, registry 1 dòng idempotent); `--help` thoát 0.

## Tiêu chí xong
- `python -B -m unittest discover -s pipeline/p08_silent_sites/tests -v` sạch.
- README ghi cách chạy, ý nghĩa cột, hai kỳ vọng viết trước, giới hạn, và ghi rõ **chưa chạy dữ liệu thật**.
- Chạy thật (bước 1 của chủ dự án): trên PC/WSL2 với D01 và `--email`; kết quả ghi `data/b11/`, tổng hợp vào `research/synthesis/` với nhãn tối đa STRONG_INFERENCE; lặp với gà/sẻ vằn trước khi nâng.
