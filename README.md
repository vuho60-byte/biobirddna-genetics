# BIRDBIODNA — Dự án "DNA gốc của chim" (Primordial Avian Genome)

> **Kho này có hai nhánh, đừng nhầm:**
> - nhánh `main` — bộ công cụ di truyền màu lông (skill BioBirdDNA Genetics v0.4.0), đã dùng được cho người nuôi chim;
> - nhánh `core-genome` — **bạn đang ở đây**, dự án nghiên cứu tìm "DNA gốc" của chào mào, đang tiến hành, nhiều kết luận còn ở mức chưa xác minh.
>
> Hai nhánh cố ý tách rời: một bên là công cụ đã chốt, một bên là nghiên cứu đang chạy. Không trộn kết luận của nhánh này sang nhánh kia.

> **Trạng thái: nghiên cứu đang tiến hành, chưa bình duyệt.** Repo công khai để minh bạch quy trình.
> Mọi khẳng định mang nhãn `CONFIRMED` / `STRONG_INFERENCE` / `UNVERIFIED`; nhiều mục còn ở mức chưa xác minh.
> Một bài do model ngoài viết đã **bị chặn** vì trích dẫn sai nguồn (xem `research/synthesis/06-review-A1.md`).
> Kết quả tính toán là số đo thật, tái lập được bằng mã trong `pipeline/`, nhưng **chưa có kiểm định thống kê so với mô hình nền**.

## Kết quả mới nhất (2026-09-08, toàn bộ bộ gene)

| Chỉ số | Ngưỡng phủ 0,8 | Ngưỡng phủ 0,5 |
|---|---|---|
| Phần tử bảo tồn (363 loài chim) ánh xạ sang chào mào | 248.461 / 2.947.588 (8,43%) | 411.275 (13,95%) |
| Tổng chiều dài trên bộ gene chào mào (đã gộp trùng lặp) | 61,50 Mb — **6,00%** | — |
| Phần vùng mã hóa của chào mào nằm trong lõi | **59,8%** (giao base-level) | — |
| Làm giàu so với rải ngẫu nhiên | **9,96 lần** | — |

### Đối chứng âm (2026-09-09)

Chạy **cùng quy trình, cùng tham số** trên bộ vùng *tiến hóa nhanh* (cùng nguồn UCSC 363 loài, chỉ khác dấu độ lệch so với mô hình trung tính):

| | Bảo tồn | Tiến hóa nhanh |
|---|---|---|
| Phần tử đưa vào | 2.947.588 | 4.054.556 |
| Ánh xạ được sang chào mào | 248.461 (8,43%) | **1.133 (0,03%)** |

Chênh **281 lần** dù bộ tiến hóa nhanh có nhiều phần tử hơn. Điều này loại bỏ giả thuyết "công cụ căn khớp bừa", nhưng **không** chứng minh vùng bảo tồn có chức năng: vùng tăng tốc khó căn vì chúng phân kỳ, mà phân kỳ là định nghĩa của chúng. Đối chứng âm đúng nghĩa (hoán vị vị trí ngẫu nhiên) chưa làm. Xem `docs/05` Phần D.

Chi tiết và 10 hướng mở rộng: `docs/05-phat-hien-va-huong-mo-rong.md`

Số liệu tính bằng `pipeline/p02_core_map/overlap_stats.py` (giao từng base, có gộp trùng lặp).
Bảng phân loại của `annotate` cộng nguyên chiều dài mỗi vùng nên cho số cao hơn (76,8% / 12,6 lần); **không dùng bảng đó cho phát biểu dạng phần trăm**.

Xếp theo tổng bp: Znf521, Sox6, Foxp1, Vps13b, Klhl29.
Xếp theo **mật độ** (bp lõi/kb locus): **Hoxa5, Hoxa11, Hoxa3**, Gpr19, Dolk, Nipbl — hơn 2.000 bp/kb.
Hai bảng khác nhau vì xếp theo tổng bp thiên vị gene dài; đọc cả hai (`pipeline/p02_core_map/README.md`).

Chi tiết số liệu: `data/a2/SUMMARY.md` · nhật ký hằng ngày: `tasks/checkpoint-*.md`

---

Nhánh nghiên cứu mở rộng từ skill **BioBirdDNA Genetics v0.4.0** (V4: hỗ trợ quyết định di truyền màu lông chim cảnh).
V4 trả lời "màu lông truyền thế nào". Dự án này hỏi câu lớn hơn: **phần nào của bộ gene chim KHÔNG đổi qua hàng chục triệu năm, phần nào đổi, và tái dựng "bộ gene tổ tiên" bằng tính toán được tới đâu.**

## Ba lớp của mô hình làm việc

| Lớp | Nội dung | Ai xử lý |
|---|---|---|
| LÕI (core-invariant) | Vùng bảo tồn cực cao qua mọi loài chim + ngoài nhóm (UCE/CNE), yếu tố điều hòa riêng lớp chim (ASHCE), gene phát triển cơ thể | Dự án này |
| VỎ (variable shell) | Màu lông, cỡ mỏ, kích thước… đột biến + chọn lọc trong vài thế hệ | V4 (`biobirddna-genetics`) |
| CÁ TÍNH (behaviour) | Hành vi bẩm sinh (di cư, hót, làm tổ) + cá tính đa gene + môi trường | Dự án này (nhánh hành vi) |

## Bản đồ thư mục

```
README.md                      ← bạn đang đọc
CLAUDE.md                      ← luật project (nhỏ, cố định)
tasks/todo.md                  ← kế hoạch theo pha, tick được
tasks/lessons.md               ← #LESSONS
tasks/primordial-dna-think.md  ← THINK (mục tiêu thật, giả định, rủi ro)
tasks/primordial-dna-plan.md   ← PLAN (tiêu chí chấp nhận, task, skill)
tasks/primordial-dna-handoff.md← bảng bàn giao subagent (append-only)
docs/00-tai-dinh-khung-cau-hoi.md  ← 6 câu hỏi gốc → câu hỏi khoa học kiểm chứng được + sổ giả thuyết H1/H2/H3
docs/01-thuat-ngu.md               ← từ điển thuật ngữ Việt–Anh
docs/02-nguon-du-lieu-va-cong-cu.md← registry dataset/tool (trạng thái: cần xác minh)
research/briefs/               ← brief giao cho subagent (R1–R6)
research/raw/                  ← kết quả thô từ subagent (không sửa tay)
research/synthesis/            ← tổng hợp của Main sau mỗi vòng
data/registry.csv              ← mọi file dữ liệu tải về: nguồn, phiên bản, sha256
pipeline/README.md             ← 6 module code dự kiến (chưa viết code khi chưa có dữ liệu)
env/environment.yml            ← môi trường conda (WSL2) cho phần tính toán
```

## Luật bằng chứng (kế thừa V4)

Mọi khẳng định trong `research/` và `docs/` mang nhãn: `CONFIRMED` (≥2 nguồn độc lập) · `STRONG_INFERENCE` (1 nguồn) · `UNVERIFIED` (chưa tìm được nguồn, nói rõ). Không bịa DOI/URL. Không nâng ngôn ngữ bằng chứng.
