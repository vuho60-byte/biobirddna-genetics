# BioBirdDNA Genetics Skill v0.3.0

Bộ khung skill và code mẫu cho dự án **Biologicalbird Community / BioBirdDNA**, tập trung vào:

- mô phỏng công thức ghép màu lông theo mô hình di truyền;
- suy luận mô hình di truyền ứng viên từ phả hệ và toàn bộ đời con;
- xử lý kiểu gene cha mẹ chưa chắc chắn bằng phân bố xác suất;
- chuẩn hóa pairing, clutch, mọi trứng/con và timeline phenotype mà không đồng nhất missing với zero;
- đo recovery, uncertainty, calibration theo design và absolute fit trên dữ liệu synthetic;
- tính quan hệ và hệ số cận huyết đời con;
- xếp hạng cặp ghép theo mục tiêu nghiên cứu, không chỉ theo màu;
- sàng lọc kiểu hình sắc tố T02 bằng rule-based code;
- giữ ranh giới rõ giữa ảnh, phenotype, mô hình di truyền và gene xác nhận.

## Thành phần chính

```text
BioBirdDNA_Genetics_Skill_v0.3.0/
├── SKILL.md
├── agents/openai.yaml
├── references/
├── planning/
├── knowledge/T02_knowledge_base_v1.0.json
├── assets/schemas/
├── assets/templates/
├── examples/input/
├── examples/output/
├── scripts/
├── src/biobirddna_genetics/
├── npm/
├── standalone/
└── tests/
```

## Tài liệu kiến trúc cho model/coder

- `BioBirdDNA_Genetics_Skill_Blueprint_v0.1.0.md`: bản hợp nhất gần 10.000 từ gồm skill contract, kiến trúc lớp, data contract, công thức, API/database draft, roadmap và handoff.
- `planning/MODEL_CODER_HANDOFF.md`: prompt, read order và quy tắc patch dành cho coding model tiếp theo.
- `planning/backlog.json`: backlog có phase, priority, acceptance criteria và dependency.
- `planning/adr/`: các quyết định kiến trúc không được thay đổi âm thầm.

## Bốn lựa chọn code

### Bản modular

Dùng để phát triển thành skill production:

```text
src/biobirddna_genetics/
```

Các module đã có:

- `cross.py`: Punnett engine cho genotype chắc chắn hoặc posterior;
- `inference.py`: likelihood/Bayes thăm dò giữa nhiều mô hình;
- `pedigree.py`: ma trận quan hệ và cận huyết;
- `optimizer.py`: xếp hạng cặp và expected information gain;
- `pigment.py`: xếp hạng phenotype T02, không xác nhận gene;
- `breeding_records.py`: chuẩn hóa hồ sơ sinh sản theo giai đoạn và đối soát số lượng;
- `simulation.py`: harness recovery/calibration synthetic có seed và cổng misspecification;
- `evidence.py`: chặn claim vượt evidence;
- `service.py`, `cli.py`: orchestration và CLI.

### Bản một file

```text
standalone/biobirddna_genetics_standalone.py
```

Bản này là snapshot tiện dụng cho `cross`, `infer` và `pigment` bằng Python
standard library. Python modular mới là implementation chuẩn cho đủ sáu operation
v0.3.0.

### Bản Node.js một file

```text
standalone/biobirddna_genetics_node.mjs
```

Bản Node.js một file vẫn là cross simulator tối giản. Package ESM trong `npm/`
không thay thế file standalone này.

### Package Node ESM

```text
npm/
```

Package không có runtime dependency và triển khai đủ sáu operation theo cùng
contract v0.3.0 với Python modular. Nó còn có feedback adapter review-only; test
packed-install xác minh cả hai CLI hoạt động ngoài cây source. Simulation dùng
cùng seed/PRNG và có cross-runtime parity ở mức parsed-JSON với float tolerance
(rel_tol `1e-9`, abs_tol `1e-12`); byte-exact chỉ được bảo đảm cho double-run
trên cùng một máy, vì libm khác nhau giữa platform có thể lệch ULP ở chữ số
thập phân cuối. Inference khóa count quan sát ở `1_000_000` để giữ work budget
và numerical parity.

## Yêu cầu

- Python 3.11 trở lên.
- Node.js 20 trở lên cho package npm/Node tùy chọn.
- Không bắt buộc thư viện ngoài.

## Cài package tùy chọn

Không cần cài để chạy các script. Khi muốn dùng entry point package:

```bash
python -m pip install --no-build-isolation -e .
biobirddna-genetics cross --input examples/input/cross_autosomal_recessive.json --pretty
```

`--no-build-isolation` hữu ích trong môi trường offline đã có sẵn setuptools.

## Package npm/Node private

Yêu cầu Node.js 20 trở lên:

```bash
cd npm
npm test
npm pack --dry-run
```

Package là `@biobirddna/genetics`, `UNLICENSED`, không có runtime dependency và
được giữ private và không cấu hình public registry scope. Đây là gói local đã được kiểm
tarball: không có npm login/publish, GitHub push hoặc xác nhận registry trong cây
skill này. Chỉ sau hành động external được owner ủy quyền mới xác minh bằng
`npm view` package/version/integrity. Hai CLI:

```text
biobirddna-genetics
biobirddna-feedback
```

Sau khi tạo tarball đã kiểm chứng, có thể cài local:

```bash
npm install --global ./biobirddna-genetics-0.3.0.tgz
biobirddna-genetics cross --input ./cross.json --output ./result.json
biobirddna-feedback --input ./review-envelope.json \
  --expected-integrity <sha256> --output ./feedback.json
```

CLI feedback chỉ tạo artifact review/proposal bất biến cho quy trình AutoResearch.
Nó không gọi experiment promotion, không tự approve/apply và không dùng health
gate 3.0 của Page làm điểm BioBirdDNA. Input là normalized production-review
wrapper, không phải raw operation result. Mutation intent, claim gene/variant
được xác nhận dưới G5 và AutoResearch health-score trá hình đều fail closed; xem
`references/autoresearch-feedback.md`.

## OAuth MCP và plugin ChatGPT/Codex cá nhân

GEN-015 thêm adapter cloud stateless trong dedicated Supabase project. Endpoint
chỉ expose đúng sáu compute tool và một feedback tool, tất cả read-only. OAuth
2.1 kiểm đúng MCP audience và allowed user trước dispatch; adapter không tạo
table, bucket, queue, model-run persistence hoặc quyền welfare override.

ChatGPT Web dùng private custom MCP App đã kết nối OAuth. Universal plugin cục bộ
cho ChatGPT desktop Work mode/Codex chỉ tham chiếu đúng registered app ID và bundle
skill đã kiểm checksum. Token, password, app ID thật và personal marketplace tree
không thuộc release source. Public authenticated review vẫn là cổng external của
GEN-016, không phải trạng thái đã hoàn tất.

Node updater/watcher chạy full source/cloud validation, phân loại sensitive diff,
yêu cầu security pin đã duyệt, deploy immutable snapshot, verify live và rollback
nếu thất bại. Watcher không gọi hoặc mutate AutoResearch.

## Kiểm chứng và cài skill global

Chạy gate đầy đủ trên đúng cây file chuẩn bị cài:

```bash
python -B scripts/validate_package.py
python -B scripts/verify_examples.py
python -B scripts/verify_package.py
```

Trên Windows, cài theo cơ chế staging, kiểm checksum rồi mới atomic rename:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/install_skill.ps1 `
  -DestinationRoot "$env:USERPROFILE\.codex\skills"
```

Trên Linux/macOS:

```bash
./scripts/install_skill.sh "${CODEX_HOME:-$HOME/.codex}/skills"
```

Installer giữ install lock, từ chối ghi đè/đua race với target đã tồn tại,
kiểm checksum trước và sau final rename, và không đưa `__pycache__`, `.pyc`
hoặc `.pyo` vào bản cài. Trên Windows có locale không mặc định UTF-8, chạy
validator metadata của Codex bằng `python -X utf8 .../quick_validate.py <skill-root>`.
Shell installer fail-loud nếu Linux/filesystem không hỗ trợ
`renameat2(RENAME_NOREPLACE)` hoặc macOS volume không hỗ trợ
`renamex_np(RENAME_EXCL)`; không fallback sang rename có thể ghi đè.

## Chạy nhanh

### Windows

```bat
run_tests.bat
run_demo.bat
```

### Linux/macOS

```bash
./run_tests.sh
./run_demo.sh
```

### Chạy từng chức năng

```bash
python -B scripts/biobirddna_cli.py cross \
  --input examples/input/cross_autosomal_recessive.json \
  --output examples/output/cross_autosomal_recessive.json

python -B scripts/biobirddna_cli.py infer \
  --input examples/input/infer_z_linked_case.json \
  --output examples/output/infer_z_linked_case.json

python -B scripts/biobirddna_cli.py pigment \
  --input examples/input/pigment_progressive_greying.json \
  --output examples/output/pigment_progressive_greying.json

python -B scripts/biobirddna_cli.py rank-mates \
  --input examples/input/mate_ranking_case.json \
  --output examples/output/mate_ranking_case.json

python -B scripts/biobirddna_cli.py normalize-records \
  --input examples/input/normalize_breeding_records.json \
  --output examples/output/normalize_breeding_records.json

python -B scripts/biobirddna_cli.py simulate-recovery \
  --input examples/input/simulate_recovery.json \
  --output examples/output/simulate_recovery.json
```

## Ví dụ công thức

### Lặn trên NST thường

```text
N/m × m/m
```

Kỳ vọng mô hình:

```text
50% N/m — phenotype bình thường nhưng có thể mang allele
50% m/m — phenotype mục tiêu nếu penetrance = 1
```

### Lặn liên kết Z ở chim

```text
Trống m/m × mái N/W
```

Kỳ vọng mô hình:

```text
50% con trống N/m — phenotype bình thường, mang allele
50% con mái m/W — phenotype mục tiêu
```

Đây chỉ là kết quả của **mô hình ký hiệu**. `N` và `m` không phải tên gene chào mào đã xác nhận.

## Kiến trúc đọc theo thứ tự

Cho model/coder mới:

1. `SKILL.md`
2. `references/scientific-boundaries.md`
3. `references/architecture.md`
4. `references/data-contracts.md`
5. `references/algorithms.md`
6. `planning/MODEL_CODER_HANDOFF.md`
7. `planning/backlog.json`
8. code trong `src/`
9. test trong `tests/`

## Phạm vi v0.3.0

Đã chạy được:

- một locus, hai allele;
- autosomal recessive/dominant/incomplete dominant;
- avian Z-linked recessive/dominant;
- penetrance đơn giản;
- lỗi chấm phenotype đơn giản;
- posterior kiểu gene cha mẹ;
- model comparison thăm dò;
- pedigree relationship matrix;
- mate ranking đa mục tiêu;
- rule-based phenotype screen.

Chưa production-ready cho:

- hai locus/epistasis;
- polygenic trait;
- survival model theo genotype;
- missing-not-at-random;
- genomic relationship matrix từ SNP;
- image segmentation/calibration trực tiếp;
- spectrophotometry/Raman/HPLC/LC-MS;
- database và API multi-user;
- clinical or invasive sampling workflow.

Các phần này đã có backlog và interface mở rộng trong `references/implementation-roadmap.md`.

## Nguyên tắc khoa học cốt lõi

```text
Ảnh → mô tả kiểu hình
Phả hệ + đời con → mô hình di truyền ứng viên
DNA đồng phân ly → locus/biến thể được hỗ trợ
Tái lập/chức năng → cơ chế xác nhận mạnh
```

Skill không cho phép nhảy thẳng từ ảnh sang gene.

## Local distribution policy (2026-09-05)

The npm package is private (`private: true`), UNLICENSED and not published. MANIFEST mirrors npm metadata; there is no public publishConfig. Only root `SYNC-NOTE.md` is a local untrusted receipt, excluded from release inventories. `verify_package.py` reports its bytes and SHA-256 separately; this digest is not authentication. Other unexpected files remain errors. Git metadata is excluded by the inventory builder.
