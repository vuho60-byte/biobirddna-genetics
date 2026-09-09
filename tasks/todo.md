# TODO — BIRDBIODNA "DNA gốc của chim"

## Trạng thái 2026-09-08 (sau sự cố xóa + khôi phục)
- [x] Khôi phục project từ transcript (66 file), 22 unittest PASS, git init + 6 commit
- [x] Đội hình mới: CẤM Codex (phải xin duyệt); Claude team + Antigravity (session `birdbiodna`)
- [x] Vá lỗi `.params` bỏ qua bước sai (nguyên nhân kết quả toàn genome hôm 07/09 giữ số liệu chr1)
- [x] A1 (Antigravity): diễn giải 10 gene lõi + 3 cảnh báo phương pháp → `research/raw/A1-*.md`
- [x] `annotate` thêm mật độ bp/kb + bảng xếp theo mật độ (+4 test) theo cảnh báo A1
- [ ] Tải lại dữ liệu lõi (~1,3 GB) — đang chạy, còn bộ gene chào mào
- [ ] V2 (Sonnet): review chéo nhãn bài A1 → `research/synthesis/06-review-A1.md`
- [x] code-x4 (Sonnet): bộ căn LAST cho phần tử <100 bp → `pipeline/p02b_last/`
- [ ] Chạy A2 toàn genome (`bash pipeline/p02_core_map/launch_full.sh`) → SUMMARY + sha256
- [ ] Báo cáo LÕI chào mào v0 theo brief S3 → `research/synthesis/03-loi-chao-mao-v0.md`
- [x] B1 2026-09-09: chạy LAST thật (23 mảnh, 13 phút), gộp BED → `data/a2_last/bulbul_core_all.bed`
      (343.781 vùng, 66,63 Mb, 6,50% bộ gene). Kết quả + giới hạn: `research/synthesis/08-ket-qua-b1-doan-ngan.md`
- [ ] Sửa bảng gene: lọc gene có CDS/locus > 0,9 (mô hình gene cụt — xem docs/05 mục D1)
- [ ] B2: kiểm định bằng hoán vị vị trí giữ nguyên tỉ lệ GC (thay cho mô hình nhị thức đã bị loại)
- [ ] Cập nhật README + docs/04 + docs/05 theo số mới của bản đồ đầy đủ
- [ ] Chạy lại sản phẩm C (20 locus, `data/c_loci` đã mất trong sự cố)

## Pha 1 — Môi trường & kiểm kê dữ liệu (chi tiết 7 việc P1.1–P1.7: research/synthesis/00-tong-hop-vong-1.md mục 9; P1.4 và P1.5 chạy song song, là GATE của Pha 2)
- [ ] WSL2 + conda env từ env/environment.yml; kiểm `python -c "import Bio"`
- [ ] Xác minh D01 (assembly Pycnonotus) — tải, BUSCO, ghi registry
- [ ] Xác minh D02/D03 (B10K alignment + track bảo tồn) — tải phần cần, ghi registry
- [ ] Quyết định: dùng track sẵn hay tự align (ghi ADR trong docs/)

## Pha 2 — Bản đồ LÕI bất biến
- [ ] Trích vùng bảo tồn ≥ ngưỡng trên ≥90% loài + ngoài nhóm → BED
- [ ] Giao với bộ gene chào mào; thống kê % genome, phân bố theo NST/vi NST
- [ ] Giao với ASHCE (D08) và UCE (D07)

## Pha 3 — Tái dựng tổ tiên (ASR) cho locus chọn
- [ ] Tập locus: gene màu lông V4 (T22) + gene hành vi Q2b
- [ ] MAFFT/PRANK → IQ-TREE2 -asr tại 3 nút (Neornithes, Passeriformes, Pycnonotidae)
- [ ] Báo cáo: vị trí khác biệt tổ tiên ↔ chào mào hiện đại

## Pha 4 — Quét chọn lọc (Q1b)
- [ ] dN/dS (codeml/HyPhy) nhóm toolkit vs nhóm hành vi vs nhóm sắc tố
- [ ] Kiểm định thống kê khác biệt giữa nhóm

## Pha 5 — Nối V4
- [ ] Đầu vào cho tracker T06 (khả thi nhánh gene) và T22 (registry gene ứng viên có nguồn/loài/mức bằng chứng)
- [ ] Không nâng ngôn ngữ bằng chứng; mọi claim gene ≤ G-level dữ liệu cho phép

## Pha 6 — Công bố & sổ giả thuyết
- [ ] Bảng H1/H2/H3 với kết quả đo thực tế cho từng dự đoán
- [ ] Bài tổng hợp tiếng Việt cho cộng đồng + bản kỹ thuật
