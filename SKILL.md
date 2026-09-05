---
name: biobirddna-genetics
description: "Use when mô phỏng phép ghép màu lông chim, chuẩn hóa hồ sơ sinh sản, đo simulation recovery/calibration, suy luận mô hình di truyền ứng viên từ phả hệ và đời con, tính quan hệ/cận huyết, xếp hạng cặp nghiên cứu hoặc sàng lọc kiểu hình sắc tố bằng code cho dự án BioBirdDNA."
---

# BioBirdDNA Genetics

Skill này là **bộ máy hỗ trợ quyết định nghiên cứu**, không phải máy xác nhận gene và không phải lệnh tự động ghép chim. Nó tách năm lớp:

1. quan sát màu nhìn thấy;
2. giả thuyết kiểu hình T02;
3. mô hình di truyền ứng viên;
4. mô phỏng đời con;
5. bằng chứng phân tử.

`<skill-root>` trong mọi lệnh là thư mục tuyệt đối chứa file `SKILL.md` này. Không diễn giải `scripts/...`, `examples/...` hoặc `knowledge/...` theo thư mục workspace khác.

## Ranh giới không được vượt

- Không chuyển tên cộng đồng như “bạch”, “bông”, “mơ”, “xám”, “khói”, “nâu” thành genotype hoặc tên gene.
- Không suy trội, lặn, liên kết Z, carrier hoặc kiểu gene từ ảnh đơn lẻ.
- Module sắc tố chỉ xếp hạng **kiểu hình nhìn thấy**; luôn giữ `biochemical_pigment_status:"unconfirmed"`, `gene_status:"not_assessed"` và `inheritance_status:"not_assessed"`.
- Không gọi xác suất model là xác suất gene đã tồn tại. `posterior` chỉ biểu thị mức phù hợp tương đối giữa các mô hình đã đưa vào so sánh dưới giả định dữ liệu và prior hiện có.
- Không gọi mô phỏng Punnett là tỷ lệ bảo đảm của một lứa nhỏ. Mỗi trứng/con là một lần lấy mẫu ngẫu nhiên, còn sống sót và chấm phenotype có thể làm lệch số quan sát.
- Không tự động đề xuất hoặc chấp thuận cặp cha–con, mẹ–con, anh chị em ruột. Mặc định chặn quan hệ trực hệ và full sibling; mọi ngưỡng cận huyết khác phải được cấu hình, ghi lý do và review.
- Không tối đa hóa màu hiếm bằng cách bỏ qua sức khỏe, phúc lợi, khả năng sinh sản, nguồn gốc, consent hoặc đa dạng di truyền.
- Không dùng skill để hướng dẫn lấy máu, sinh thiết, gây mê hoặc thao tác xâm lấn. Mẫu phân tử thuộc SOP thú y/phòng xét nghiệm riêng.
- Không ghi `confirmed_gene`, `causal_variant` hoặc “cơ chế đã xác nhận” dưới mức G5; G6 mới cho phép ngôn ngữ xác nhận mạnh.
- Không sửa dữ liệu gốc để làm khớp mô hình. Bản ghi đời con phải gồm cả bình thường, bất thường, chưa xác định, chết phôi, không nở và mất theo dõi khi có dữ liệu.
- Không loại `M0_NONGENETIC` hoặc `UNK` chỉ để hệ thống luôn trả một gene.
- Không dùng nhiều ảnh cùng một cá thể như nhiều cá thể độc lập khi huấn luyện/đánh giá model.

## Mặc định an toàn

- Invocation ngầm chỉ đọc dữ liệu cục bộ, kiểm định, mô phỏng và xuất report. Không ghi vào Drive, Supabase, Google Sheet, phòng lab hoặc hệ thống bên ngoài.
- Mọi kết quả ghép có trạng thái `decision_support_only` cho đến khi reviewer người duyệt.
- Founder/cha mẹ không có trong pedigree là **quan hệ chưa biết**, không phải bằng chứng không liên quan; `offspring_inbreeding:0.0` của prototype không được diễn giải là “đã chứng minh không cận huyết”.
- `health_risk` hiện là thành phần trừ điểm; chỉ `welfare_blocked` và các block phả hệ/cận huyết được mã hóa thành hard gate. Không dùng score để xác nhận một cặp là an toàn sức khỏe.
- Thiếu giới tính đáng tin cậy thì không dùng kết quả để phân biệt mô hình liên kết Z.
- Thiếu số con đầy đủ thì hạ evidence grade; không tự điền số con bị thiếu bằng 0.
- Thiếu T02 phenotype hoặc ảnh đạt Q2 trở lên thì module sắc tố trả `UNK`/`needs_recapture` thay vì đoán.

## Định tuyến yêu cầu

### 1. Mô phỏng phép ghép đã biết mô hình

Dùng khi người dùng cung cấp mô hình và kiểu gene cha mẹ hoặc phân bố xác suất kiểu gene.

Lệnh:

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" cross --input <absolute-cross-json> --output <absolute-output-json>
```

Đầu ra bắt buộc:

- mô hình đang giả định;
- kiểu gene hoặc posterior cha mẹ;
- xác suất kiểu gene đời con;
- xác suất phenotype theo giới;
- giới hạn `simulation_not_confirmation`.

### 2. Suy luận mô hình di truyền ứng viên

Dùng khi có phenotype bố mẹ và số đời con theo giới/phenotype.

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" infer --input <absolute-inference-json> --output <absolute-output-json>
```

Luôn so sánh nhiều mô hình, tối thiểu:

- `M1_AR` — lặn trên NST thường;
- `M2_AD` — trội trên NST thường;
- `M3_INC` — trội không hoàn toàn;
- `M4_ZR` — lặn liên kết Z;
- `M5_ZD` — trội liên kết Z;
- `M0_NONGENETIC` — nền không di truyền.

Nếu caller truyền candidate set tùy chỉnh, runtime yêu cầu ít nhất hai model và
đúng một `M0_NONGENETIC`; trọng số prior là số hữu hạn, không âm và chỉ có nghĩa
tương đối. `model_id` phải duy nhất. Không dùng một candidate duy nhất để tạo
posterior `1.0`.

Mọi inheritance model phải có `normal_label != target_label`; riêng
`incomplete_dominant`, `intermediate_label` phải khác cả hai label còn lại.
Mỗi count offspring là integer trong `[0, 1_000_000]`.

Không báo “gene lặn” chỉ vì mô hình lặn đứng hạng nhất. Câu đúng là:

> Dữ liệu hiện tại phù hợp hơn với mô hình … trong tập mô hình đã kiểm tra; cần các phép ghép/đời con/DNA bổ sung để loại trừ mô hình cạnh tranh.

### 3. Xếp hạng cặp ghép nghiên cứu

Dùng khi có:

- danh sách chim trống/mái;
- posterior kiểu gene theo từng mô hình;
- pedigree;
- model weights;
- risk flags.

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" rank-mates --input <absolute-ranking-json> --output <absolute-output-json>
```

Thiết kế production phải tách hai mục tiêu:

- `research_pairing`: cặp giúp phân biệt mô hình tốt nhất;
- `target_pairing`: cặp tăng xác suất phenotype mục tiêu nhưng vẫn qua cổng phả hệ, sức khỏe và phúc lợi.

Prototype v0.3.0 hiện trả một ranking có trọng số chung và chưa triển khai tách
hai objective này; đó là GEN-040 trong backlog. Vì vậy không được gắn nhãn
ranking hiện tại là `research_pairing` hay `target_pairing`, và không dùng nó
như phê duyệt ghép thật.

Hypothesis weights và pairing weights là trọng số tương đối hữu hạn, không âm;
hypothesis weights được normalize. Runtime fail closed nếu pair score tràn thành
số không hữu hạn.

Một cặp bị block vẫn được giữ trong report với lý do; không xóa để làm bảng đẹp hơn.

### 4. Sàng lọc kiểu hình sắc tố

Dùng khi có feature đã đo từ ảnh/quan sát T02, không dùng trực tiếp file ảnh ở v0.3.0.

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" pigment --input <absolute-pigment-json> --output <absolute-output-json>
```

Đầu ra chỉ dùng các nhãn làm việc:

`MEL-LEU`, `MEL-PGR`, `MEL-ALB`, `MEL-INO`, `MEL-BRN`, `MEL-DIL`, `MEL-MEL`, `CAR-ALT`, `LOC-ACQ`, `IMG-ART`, `BASE-NORM`, `UNK`.

Không dùng đầu ra này để tạo công thức ghép cho đến khi có dữ liệu phả hệ/đời con độc lập.

### 5. Chuẩn hóa hồ sơ sinh sản

Dùng khi cần đưa pairing, clutch, số trứng theo giai đoạn và từng hồ sơ
trứng/con về một contract có thể audit:

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" normalize-records --input <absolute-record-json> --output <absolute-output-json>
```

Bắt buộc giữ:

- nguồn xác định giới tính và confidence của cả hai bố mẹ;
- `missing` khác `observed zero`;
- toàn bộ trứng/con, gồm bình thường, không thụ tinh, chết phôi, không nở,
  chết sau nở và mất theo dõi;
- chiều trống × mái của reciprocal cross;
- đối soát stage counts khi coverage của hồ sơ con là `complete`;
- phenotype timeline theo `sequence_index`, không suy gene từ timeline.

### 6. Đo recovery trên dữ liệu synthetic

Dùng để kiểm tra một candidate-model set có khôi phục được mô hình sinh dữ liệu
đã biết hay không:

```text
python -B "<skill-root>/scripts/biobirddna_cli.py" simulate-recovery --input <absolute-simulation-json> --output <absolute-output-json>
```

Harness bắt buộc:

- giữ đúng một `M0_NONGENETIC` và ít nhất một model đối chứng khác;
- dùng genotype đúng `NA` cho generating model `non_genetic`;
- dùng seed SHA-256 53-bit và LCG chung cho Python/Node, nằm trong JSON
  safe-integer; cùng request v0.3.0 phải tạo cùng stream/output;
- giữ nguyên thứ tự `candidate_model_ids` trong output, nhưng tie set được sắp
  theo model ID thay vì dùng thứ tự input để phá tie;
- báo unique top-1 và tie set tách biệt;
- chỉ gọi calibration khi design có nhiều realized true models;
- refit parent configuration trong từng parametric-bootstrap replicate;
- cảnh báo misspecified set từ absolute fit, không từ posterior tương đối;
- trả `synthetic_only`, `decision_support_only`, `gene_status:not_assessed`;
- không mô hình hóa genotype-specific survival trước GEN-012.

Built-in model IDs giữ semantics canonical ngay cả khi generating model không
nằm trong candidate set. Generator misspec tùy biến phải dùng ID riêng.

### 7. Dùng package Node và tạo feedback review-only

Package ESM trong `npm/` cung cấp API/CLI parity cho sáu operation trên Node.js
>=20, không cần Python và không có dependency runtime:

```text
biobirddna-genetics <operation> --input <absolute-json> --output <absolute-json>
```

Để chuyển một production envelope đã hoàn tất thành artifact feedback:

```text
biobirddna-feedback --input <absolute-envelope-json> --expected-integrity <sha256> --output <absolute-artifact-json>
```

Feedback adapter chỉ tạo artifact bất biến `pending_human_review`,
`auto_apply:false`, `target_files:[]`. Nó không gọi AutoResearch, không dùng
health gate 3.0 làm score BioBirdDNA, không ghi proposal ledger và không tự apply
thay đổi. Mutation intent, claim gene/variant đã xác nhận dưới G5 và
AutoResearch health-score trá hình đều bị chặn kể cả camel/Pascal case. Đọc
contract đầy đủ tại `references/autoresearch-feedback.md`.

### 8. Dùng OAuth MCP trên ChatGPT Web hoặc plugin cá nhân

Remote MCP của release này chỉ là transport stateless quanh bảy tool cố định.
OAuth phải kiểm exact resource audience và đúng allowed user trước khi tool chạy.
Không có anonymous fallback, generic dispatch, database result store, RBAC cho
stored records, welfare override hoặc publication flow.

Sáu tool đầu gọi đúng sáu operation Node hiện có. Tool feedback thứ bảy chỉ tạo
artifact review-only với integrity do server sở hữu; output vẫn
`pending_human_review`, `auto_apply:false`, `target_files:[]`. Plugin cá nhân chỉ
tham chiếu connection ChatGPT đã đăng ký và bundle cây skill đã kiểm checksum;
không bundle secret/token hoặc thêm MCP server thứ hai.

### 9. Audit dữ liệu hoặc code

- Đọc `references/scientific-boundaries.md` trước.
- Đọc `references/data-contracts.md` để kiểm schema.
- Chạy test:

```text
python -B -m unittest discover -s "<skill-root>/tests" -v
```

- Chạy kiểm toàn gói:

```text
python -B "<skill-root>/scripts/validate_package.py"
```

Không công nhận patch nếu test cũ bị xóa, bỏ qua hoặc đổi expected result chỉ để pass.

## Báo cáo sau khi chạy CLI

CLI modular là primitive tính toán của prototype: envelope hiện tại chứa
`ok`, `operation`, `engine_version`, `result` và `result.limits`. Nó không tự
chứng nhận evidence grade và không thay thế run envelope production trong
`references/data-contracts.md`.

Mọi câu trả lời cho người dùng sau khi chạy CLI phải bổ sung hai mục rõ ràng:

- **Evidence:** nguồn input đã dùng; `G0`–`G6` chỉ khi hồ sơ cung cấp đủ căn cứ
  theo bảng evidence bên dưới, nếu không ghi `not_assessed`; posterior cao
  không tự nâng evidence grade.
- **Limits:** lặp lại giới hạn từ `result.limits`, dữ liệu thiếu, trạng thái
  `decision_support_only` và yêu cầu human review khi có pairing.

Không tự bịa `run_id`, `input_sha256`, evidence grade hoặc
`human_review_status` nếu pipeline production chưa cung cấp các trường đó.

## Cổng dữ liệu trước khi chạy

### Gate P — Phenotype

Phải có tối thiểu:

- `individual_id` ổn định;
- tên cộng đồng giữ nguyên ở nguồn;
- phenotype chuẩn hóa theo T02;
- tuổi/giai đoạn lông;
- Q-level;
- cùng cá thể qua thời gian nếu kết luận phụ thuộc tiến triển.

### Gate G — Genetics

Để chạy inference có ý nghĩa cần:

- giới tính từng parent và tốt nhất từng con;
- phương pháp xác định giới;
- ID cặp ghép/lứa;
- tổng trứng, phôi, nở, sống và phenotype theo mốc;
- cả con bình thường và con màu;
- quan hệ pedigree;
- không trộn nhiều cặp vào một dòng tổng hợp.

### Gate M — Molecular

Để nâng lên G5 cần:

- sample ID nối đúng individual ID;
- chain of custody;
- variant/locus;
- genotype từng cá thể;
- cosegregation phenotype–genotype;
- QC phòng xét nghiệm;
- phân tích độc lập với ảnh nếu có thể.

## Trình tự chuẩn

1. Đọc `knowledge/T02_knowledge_base_v1.0.json` và giữ nguyên ontology T02.
2. Xác định yêu cầu là `cross`, `infer`, `rank-mates`, `pigment` hay `audit`.
3. Kiểm cổng dữ liệu; thiếu gì phải báo trước khi suy luận.
4. Chạy deterministic core trước. Không gọi AI/ML cho việc Punnett có thể tính chính xác.
5. Với inference, marginalize genotype cha mẹ chưa biết; không ép một cấu hình duy nhất.
6. Với ranking, tính cận huyết và block relation trước khi tính điểm hấp dẫn màu.
7. Luôn xuất bất định, giả thuyết thay thế và giới hạn.
8. Gắn G-level cho claim di truyền.
9. Reviewer người duyệt trước mọi quyết định ghép thật.
10. Khi có lứa mới, append dữ liệu rồi chạy lại; không sửa số cũ.

## Hợp đồng kiểu gene

### Locus thường

```text
N/N
N/m
m/m
```

### Locus liên kết Z ở chim

```text
Trống ZZ: N/N, N/m, m/m
Mái ZW:   N/W, m/W
```

`N` và `m` chỉ là allele ký hiệu trong mô hình. Chúng không phải tên gene thật.

## Hợp đồng evidence G0–G6

| Mức | Dữ liệu | Câu chữ tối đa |
| --- | --- | --- |
| G0 | Ảnh/kiểu hình | Chưa tạo kết luận di truyền |
| G1 | Một cặp hoặc một lứa | Mô hình thăm dò |
| G2 | Nhiều lứa cùng cặp | Mô hình ứng viên |
| G3 | Nhiều cặp độc lập, có lai thuận–nghịch | Mô hình phân ly được hỗ trợ |
| G4 | Phả hệ nhiều thế hệ/test cross | Mô hình di truyền mạnh, chưa xác nhận locus |
| G5 | DNA đồng phân ly | Locus/biến thể được dữ liệu hỗ trợ |
| G6 | Tái lập độc lập hoặc bằng chứng chức năng | Cơ chế được xác nhận mạnh |

## State report hợp lệ

- `input_invalid`
- `phenotype_only`
- `insufficient_evidence`
- `simulation_ready`
- `model_inference_ready`
- `mate_ranking_ready`
- `blocked_by_pedigree`
- `blocked_by_welfare`
- `human_review_required`
- `molecular_review_required`
- `validated`
- `verification_failed`

Không gọi `simulation_ready` là `gene_confirmed`.

## Hợp đồng đầu ra tối thiểu

```yaml
operation: cross | infer | pigment | rank-mates | normalize-records | simulate-recovery
engine_version: 0.3.0
case_id: optional
pairing_id: optional
observations: {}
model_assumptions: []
results: {}
uncertainty: {}
evidence_grade: G0-G6
blocked_reasons: []
limits:
  phenotype_only: true|false
  biochemical_pigment_status: unconfirmed|measured
  gene_status: not_assessed|candidate|cosegregation_supported|replicated
  decision_support_only: true
missing_data: []
human_review_status: pending|reviewed
```

## Bản đồ tài nguyên

- Kiến trúc sâu và dependency rules: [references/architecture.md](references/architecture.md)
- Công thức và thuật toán: [references/algorithms.md](references/algorithms.md)
- JSON contract và ID: [references/data-contracts.md](references/data-contracts.md)
- Ranh giới khoa học/an toàn: [references/scientific-boundaries.md](references/scientific-boundaries.md)
- Tích hợp ontology T02: [references/t02-integration.md](references/t02-integration.md)
- Kế hoạch phát triển production: [references/implementation-roadmap.md](references/implementation-roadmap.md)
- Hướng dẫn model khác tiếp tục code: [planning/MODEL_CODER_HANDOFF.md](planning/MODEL_CODER_HANDOFF.md)
- Backlog máy đọc: [planning/backlog.json](planning/backlog.json)
- Bộ mô hình mặc định: [assets/models/default-models.json](assets/models/default-models.json)
- JSON Schema: [assets/schemas](assets/schemas)
- Request mẫu: [assets/templates](assets/templates)
- Dữ liệu demo: [examples/input](examples/input)
- Code modular: [src/biobirddna_genetics](src/biobirddna_genetics)
- Code Python một file: [standalone/biobirddna_genetics_standalone.py](standalone/biobirddna_genetics_standalone.py)
- Code Node.js một file: [standalone/biobirddna_genetics_node.mjs](standalone/biobirddna_genetics_node.mjs)
- Package Node ESM: [npm/package.json](npm/package.json)
- Contract feedback review-only: [references/autoresearch-feedback.md](references/autoresearch-feedback.md)

## Tiêu chuẩn hoàn thành một thay đổi code

Một model khác chỉ được báo hoàn thành khi:

1. nêu rõ issue/backlog ID;
2. không phá các invariant trong skill;
3. thêm hoặc cập nhật test có ý nghĩa;
4. chạy toàn bộ test sạch;
5. chạy demo liên quan;
6. cập nhật schema/docs nếu contract đổi;
7. tăng version nếu thay semantics;
8. ghi changelog;
9. không nâng ngôn ngữ bằng chứng;
10. cung cấp diff hoặc danh sách file thay đổi.
