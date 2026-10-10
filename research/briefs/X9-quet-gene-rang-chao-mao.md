# X9 — Việc B5: quét mảnh vỡ gene răng (và GULO) trên bộ gene chào mào

Chỉ tạo/sửa trong `pipeline/p07_tooth_genes/**` (+ 1 dòng `env/environment.yml`, `tasks/todo.md`, `tasks/primordial-dna-handoff.md`). **Không có lệnh xóa nào.** Python 3.11 **stdlib**. Không chạy docker/mạng ở bước viết code; BLAST+ thật chỉ dùng để chạy test tích hợp trên dữ liệu tổng hợp nhúng.

Người thực hiện: Main (phiên 2026-10-10), không giao subagent. Brief này ghi lại phạm vi để lần chạy thật sau theo đúng khung.

## Câu hỏi khoa học (docs/00 Q4; tổng hợp `research/synthesis/11`, mục 7 việc 1–2)

Nếu chim hiện đại chia sẻ **một** tổ tiên đã mất răng khoáng hóa (~116 Ma — Meredith 2014, R12 dòng #23, CONFIRMED), thì 6 gene men/ngà răng **ENAM, AMBN, AMELX, AMTN, DSPP, MMP20** phải còn trong bộ gene chào mào dưới dạng **giả gene** mang tổn thương bất hoạt (codon dừng sớm, lệch khung) ở vị trí **khớp** với các bộ chim khác. H3 (bộ gene được thiết kế riêng) không dự đoán một bộ gene "sạch" phải mang gene hỏng vô dụng cùng kiểu lỗi với họ hàng. Đây là phép thử **phân biệt giả thuyết** đầu tiên chạy được trên chính loài của dự án (docs/05 B5), và chưa ai làm ở bộ Sẻ (tổng hợp 07 mục 5).

**GULO** đi kèm: lấp khoảng trống #19 của R15 (chim chỉ có dữ liệu enzyme 1969; trạng thái gene ở sẻ chưa đọc từ trình tự).

Logic phép thử:
- Tìm thấy mảnh + tổn thương **khớp** Meredith 2014 → ủng hộ H1/H2 trực tiếp trên chào mào.
- Tìm thấy mảnh, tổn thương **không khớp** / riêng chào mào → vẫn tương thích H1 (giả gene tích thêm đột biến sau khi hỏng), không phải bằng chứng chống H1.
- **Không tìm thấy** → **không kết luận** (lắp ráp mức scaffold N50 ~218 kb; giả gene 116 Ma phân kỳ mạnh; ngưỡng e-value).
- Đối chứng dương (MEPE, SPARC — gene còn chức năng ở chim, MEPE cùng họ SCPP với gene răng) **không đạt** → toàn bộ kết quả gene mục tiêu không diễn giải được.

## Dữ liệu
- `data/ncbi/GCA_013400435.1/..._genomic.fna.gz` — bộ gene chào mào (D01, 325 MB nén). Cần giải nén (~1,2 GB) cho `makeblastdb`; **không commit**.
- Protein truy vấn ngoài nhóm: lấy từ NCBI Protein (RefSeq) bằng E-utilities, loài mặc định cá sấu Mỹ *Alligator mississippiensis* (D06 docs/02). **Không ghi accession từ trí nhớ** — `fetch-queries` tra esearch rồi ghi accession thật vào `queries_manifest.tsv` và 1 dòng `data/registry.csv` (sha256).

## Việc — `pipeline/p07_tooth_genes/tooth_scan.py`, năm lệnh con

1. `fetch-queries --genes ENAM,AMBN,AMELX,AMTN,DSPP,MMP20,GULO --controls MEPE,SPARC --organism "Alligator mississippiensis" --email X --out queries.faa --manifest queries_manifest.tsv [--registry data/registry.csv]`
   esearch `GENE[Gene Name] AND "loài"[Organism] AND refseq[Filter]` → efetch FASTA → giữ **isoform dài nhất**; header `GENE|ROLE|ACC`. Gene không có ở ngoài nhóm → dòng manifest `not_found` (không crash). Registry: thêm đúng 1 dòng theo header hiện có, idempotent, không xóa/ghi đè.
2. `parse --tblastn tblastn.tsv [--manifest ...] --out-loci loci.tsv --out-lesions lesions.tsv [--out-hits hits.tsv]`
   Đọc bảng tblastn `-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sframe qseq sseq"`; nối HSP cùng (truy vấn, scaffold, mạch) cách ≤ 100 kb thành locus; mỗi locus: độ phủ truy vấn (hợp các đoạn), e-value tốt nhất, %id, **codon dừng** (`*` trong `sseq`) chia *thô* / *có chỗ dựa* (vị trí thuộc HSP "chủ sở hữu" = bitscore cao nhất tại vị trí đó, cánh ≥10 aa mỗi bên khớp ≥50 %), **lệch khung ứng viên** (2 HSP kề nhau trên truy vấn, khác khung, khoảng cách subject ≤ 40 nt), **ranh exon** (khoảng cách subject > 40 nt — đổi khung giữa exon là bình thường).
3. `extract --fna genome.fna.gz --loci loci.tsv --pad 2000 --out loci.fa --translate loci_translate.txt` — trích vùng (theo hướng đọc gene) + 3 khung dịch để **kiểm thủ công** (yêu cầu của tổng hợp 07 mục 5).
4. `translate --fasta X --out Y` — 3 khung dịch cho FASTA bất kỳ.
5. `report --loci --lesions [--manifest] --out report.md` — bảng đối chứng (ĐẠT/KHÔNG ĐẠT), bảng gene mục tiêu, bảng tổn thương `supported=1` với cột **"Khớp Meredith 2014?"** để trống cho người điền, và **phần cảnh báo bắt buộc** (thứ tự khả năng giải thích: lỗi lắp ráp → giả tạo căn trình tự → ranh exon → tổn thương thật; "không tìm thấy ≠ không tồn tại"; nhãn tối đa STRONG_INFERENCE; cấm "đã chứng minh").

Wrapper `run_tooth_scan.sh` (Git Bash + docker như `run_last.sh`, hoặc `--mode local` cho WSL2/Linux): fetch → gunzip → makeblastdb → tblastn (`-seg no`, `-max_target_seqs 50`) → parse → extract → report → SHA256SUMS; resume theo `.params` (hàm sao chép nguyên văn từ `run_a2.sh`); `PROJECT_DIR` lên **2 cấp** (lesson 2026-09-09). Image BLAST mặc định `quay.io/biocontainers/blast:2.16.0--h66d330f_5` — tag suy từ bioconda repodata, **chưa pull**, chủ dự án duyệt trước (docs/02: "cần image BLAST, chưa duyệt").

## Test (`tests/test_tooth_scan.py`, unittest, dữ liệu nhúng)
- Bảng HSP nhúng là **đầu ra thật của tblastn 2.12.0+** trên bộ gene tổng hợp (seed cố định): nguyên vẹn / codon dừng + lệch khung mạch + / cùng locus đảo ngược mạch − / 2 exon + intron 400 nt.
- Kiểm: dịch mã 3 khung; đọc bảng + lọc e-value + sai số cột; nối locus; codon dừng đúng vị trí aa 91 và tọa độ nt 571 (+) / 567 (−); `*` giả trong HSP yếu mở rộng sai khung **không** được tính; lệch khung 1→3 / −1→−3; ranh exon **không** bị gọi là tổn thương dù HSP chồng lan 50 aa trên truy vấn; `extract` cắt pad ở biên scaffold, đảo mạch −, dịch khung +1 ra đúng protein; `report` khóa diễn giải khi đối chứng không đạt / không có đối chứng, in đủ cảnh báo; `fetch-queries` với mạng giả lập (chọn isoform dài nhất, `not_found`, registry thêm 1 dòng idempotent); `--help` mọi lệnh thoát 0; **test tích hợp** chạy lại makeblastdb + tblastn thật nếu có trên PATH.

## Tiêu chí xong
- `python -B -m unittest discover -s pipeline/p07_tooth_genes/tests -v` sạch; `bash -n run_tooth_scan.sh` sạch; smoke test `--mode local` trên bộ gene tổng hợp chạy hết 8 bước và lần 2 toàn SKIP; một lần chạy **không** truyền `--data-dir` phải báo thiếu bộ gene đúng đường dẫn mặc định (exit 4).
- README ghi: cách chạy, ý nghĩa cột, **thứ tự khả năng giải thích**, giới hạn (không phát hiện tổn thương ở vùng không căn được; lệch khung trùng ranh exon không phân biệt được; `*` do lỗi consensus lắp ráp; chưa chạy dữ liệu thật).
- Thêm 1 dòng `tasks/primordial-dna-handoff.md`, 1 dòng `tasks/todo.md`, `blast` vào `env/environment.yml`.

## Trả về ≤10 dòng
file tạo · dòng tổng kết unittest · xác nhận không lệnh xóa · giả định/giới hạn.
