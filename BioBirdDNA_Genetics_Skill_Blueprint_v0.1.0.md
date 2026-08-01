# BioBirdDNA Genetics Skill — Blueprint hợp nhất v0.1.0

Tài liệu hợp nhất dành cho model/coder tiếp tục triển khai skill. Nguồn từng phần vẫn được giữ riêng trong package để dễ version và review.

Không dùng tài liệu này để xác nhận gene từ ảnh hoặc tự động ra lệnh ghép chim.


---

<!-- BEGIN SKILL CONTRACT: SKILL.md -->

---
name: biobirddna-genetics
description: "Use when mô phỏng phép ghép màu lông chim, suy luận mô hình di truyền ứng viên từ phả hệ và đời con, tính quan hệ/cận huyết, xếp hạng cặp nghiên cứu hoặc sàng lọc kiểu hình sắc tố bằng code cho dự án BioBirdDNA."
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
- Thiếu giới tính đáng tin cậy thì không dùng kết quả để phân biệt mô hình liên kết Z.
- Thiếu số con đầy đủ thì hạ evidence grade; không tự điền số con bị thiếu bằng 0.
- Thiếu T02 phenotype hoặc ảnh đạt Q2 trở lên thì module sắc tố trả `UNK`/`needs_recapture` thay vì đoán.

## Định tuyến yêu cầu

### 1. Mô phỏng phép ghép đã biết mô hình

Dùng khi người dùng cung cấp mô hình và kiểu gene cha mẹ hoặc phân bố xác suất kiểu gene.

Lệnh:

```text
python "<skill-root>/scripts/biobirddna_cli.py" cross --input <absolute-cross-json> --output <absolute-output-json>
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
python "<skill-root>/scripts/biobirddna_cli.py" infer --input <absolute-inference-json> --output <absolute-output-json>
```

Luôn so sánh nhiều mô hình, tối thiểu:

- `M1_AR` — lặn trên NST thường;
- `M2_AD` — trội trên NST thường;
- `M3_INC` — trội không hoàn toàn;
- `M4_ZR` — lặn liên kết Z;
- `M5_ZD` — trội liên kết Z;
- `M0_NONGENETIC` — nền không di truyền.

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
python "<skill-root>/scripts/biobirddna_cli.py" rank-mates --input <absolute-ranking-json> --output <absolute-output-json>
```

Skill phải tách hai mục tiêu:

- `research_pairing`: cặp giúp phân biệt mô hình tốt nhất;
- `target_pairing`: cặp tăng xác suất phenotype mục tiêu nhưng vẫn qua cổng phả hệ, sức khỏe và phúc lợi.

Một cặp bị block vẫn được giữ trong report với lý do; không xóa để làm bảng đẹp hơn.

### 4. Sàng lọc kiểu hình sắc tố

Dùng khi có feature đã đo từ ảnh/quan sát T02, không dùng trực tiếp file ảnh ở v0.1.0.

```text
python "<skill-root>/scripts/biobirddna_cli.py" pigment --input <absolute-pigment-json> --output <absolute-output-json>
```

Đầu ra chỉ dùng các nhãn làm việc:

`MEL-LEU`, `MEL-PGR`, `MEL-ALB`, `MEL-INO`, `MEL-BRN`, `MEL-DIL`, `MEL-MEL`, `CAR-ALT`, `LOC-ACQ`, `IMG-ART`, `BASE-NORM`, `UNK`.

Không dùng đầu ra này để tạo công thức ghép cho đến khi có dữ liệu phả hệ/đời con độc lập.

### 5. Audit dữ liệu hoặc code

- Đọc `references/scientific-boundaries.md` trước.
- Đọc `references/data-contracts.md` để kiểm schema.
- Chạy test:

```text
python -m unittest discover -s "<skill-root>/tests" -v
```

- Chạy kiểm toàn gói:

```text
python "<skill-root>/scripts/validate_package.py"
```

Không công nhận patch nếu test cũ bị xóa, bỏ qua hoặc đổi expected result chỉ để pass.

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
operation: cross | infer | pigment | rank-mates
engine_version: 0.1.0
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


<!-- END SKILL CONTRACT -->

---

<!-- BEGIN ARCHITECTURE: references/architecture.md -->

# Kiến trúc sâu — BioBirdDNA Genetics Skill

**Mã kiến trúc:** `CBIM-DNA-GEN-ARCH`  
**Phiên bản:** `0.1.0`  
**Đối tượng:** model/coder tiếp tục xây skill, reviewer khoa học, người quản trị dữ liệu và người vận hành dự án.

---

## 1. Mục tiêu kiến trúc

Hệ thống phải trả lời bốn nhóm câu hỏi nhưng không trộn chúng:

1. **Màu nhìn thấy là gì?**
2. **Dữ liệu đời con phù hợp hơn với mô hình di truyền nào?**
3. **Nếu giả định một mô hình và kiểu gene cha mẹ, đời con có phân bố nào?**
4. **Cặp nào có giá trị nghiên cứu cao, ít quan hệ huyết thống hơn và vẫn qua cổng sức khỏe/phúc lợi?**

Kiến trúc phải hỗ trợ mở rộng từ prototype cục bộ sang dịch vụ nhiều người dùng nhưng giữ lõi deterministic và evidence boundary không đổi.

### 1.1. Mục tiêu chức năng

- Mô phỏng phép ghép một locus trên NST thường và Z.
- Mô phỏng cha mẹ có posterior genotype, không ép nhãn cứng.
- So sánh nhiều mô hình cạnh tranh bằng likelihood/Bayes.
- Tính ma trận quan hệ, kinship và inbreeding của đời con.
- Tính expected information gain của phép ghép nghiên cứu.
- Xếp hạng cặp theo hàm mục tiêu có hard constraints.
- Nhận feature phenotype T02 và xếp hạng nhãn làm việc.
- Xuất JSON có version, assumptions, evidence và limits.
- Cho phép model khác mở rộng qua interface ổn định.

### 1.2. Mục tiêu phi chức năng

- **Deterministic-first:** Punnett, pedigree và validation phải cho kết quả tái lập.
- **Explainable:** mọi score có thành phần và lý do.
- **Fail-closed:** thiếu dữ liệu hoặc ngoài phân bố thì trả `UNK`/blocked.
- **Versioned:** mọi model/schema/output đều có version.
- **Auditable:** giữ input hash, engine version, config, result và reviewer.
- **No hidden promotion:** không tự nâng phenotype thành genotype/gene.
- **Portable:** lõi v0.1.0 chạy bằng Python standard library.
- **Testable:** mọi invariant có unit test hoặc integration test.

### 1.3. Không thuộc kiến trúc hiện tại

- Điều khiển sinh sản tự động.
- Chẩn đoán bệnh qua ảnh.
- Xác nhận thành phần pigment hóa học từ RGB.
- Gọi gene causal chỉ từ association yếu.
- Hướng dẫn thao tác mẫu xâm lấn.
- Định giá thương mại màu hiếm.

---

## 2. Bối cảnh hệ thống

```text
Người nuôi / nghệ nhân / reviewer
            │
            ▼
     BioBirdDNA intake
            │
       ┌────┴────┐
       │         │
       ▼         ▼
 T02 phenotype  Pedigree + brood records
       │         │
       └────┬────┘
            ▼
 BioBirdDNA Genetics Skill
 ┌──────────┼──────────┬───────────┐
 ▼          ▼          ▼           ▼
Cross     Inference   Pedigree    Pigment
Engine    Engine      /Mate       Screen
                       Optimizer
            │
            ▼
  Evidence-bounded report
            │
            ▼
Human scientific / welfare review
```

### 2.1. Tác nhân

- `contributor`: gửi ảnh, phả hệ, lứa ghép.
- `data_steward`: chuẩn hóa ID và kiểm dữ liệu thiếu.
- `phenotype_reviewer`: chấm T02.
- `genetics_reviewer`: duyệt mô hình và phép ghép kiểm chứng.
- `veterinary_reviewer`: duyệt sức khỏe/phúc lợi.
- `lab_reviewer`: duyệt dữ liệu phân tử.
- `model_operator`: chạy skill, không tự thay dữ liệu nguồn.
- `publication_reviewer`: quyết định câu chữ công khai.

### 2.2. Trust boundary

- Ảnh, lời kể và seller claim là input không tin cậy mặc định.
- JSON sau validation chỉ “đúng cấu trúc”, chưa đồng nghĩa “đúng sinh học”.
- Output của code là artifact có thể audit, chưa phải quyết định cuối.
- Chỉ dữ liệu lab có chain of custody mới được dùng cho G5.

---

## 3. Cấu trúc container/component

### 3.1. Container hiện tại

```text
CLI local
  ├── JSON file adapter
  ├── deterministic genetics core
  ├── probabilistic inference
  ├── pedigree / mate optimizer
  ├── phenotype rule engine
  └── JSON report writer
```

### 3.2. Container production tương lai

```text
Web/API Gateway
    │
    ├── Auth / consent / role checks
    ├── Case Service
    ├── Genetics Compute Service
    ├── Media Feature Service
    ├── Pedigree Service
    ├── Model Registry
    ├── Audit/Event Store
    └── Report/Publication Service
```

Không đưa database, web framework hoặc queue vào deterministic core. Core chỉ nhận domain object và trả domain result.

---

## 4. Kiến trúc phân lớp

## L0 — Contracts và versioning

Chứa:

- JSON Schema;
- ID conventions;
- enum;
- model version;
- output contract;
- evidence language.

**Invariant:** thay nghĩa field phải tăng version; không reuse field cũ với semantics mới.

## L1 — Domain model

Các object chính:

- `InheritanceModel`
- `Genotype`
- `GenotypePosterior`
- `OffspringObservation`
- `PedigreeIndividual`
- `GeneticProfile`
- `ModelHypothesis`
- `PigmentObservation`
- `EvidenceAssessment`

L1 không được import adapter, database hoặc UI.

## L2 — Deterministic genetics core

Chứa:

- parser/canonical genotype;
- gamete generator;
- phenotype mapping;
- cross known genotypes;
- cross uncertain genotypes.

L2 phải là pure function nếu có thể.

## L3 — Statistical inference

Chứa:

- likelihood;
- parent genotype marginalization;
- model prior/posterior;
- observation error;
- future: survival/missingness, multi-locus, hierarchical model.

L3 được import L1–L2, không import optimizer/UI.

## L4 — Pedigree và population control

Chứa:

- topological pedigree validation;
- numerator relationship matrix;
- ancestor/sibling detection;
- offspring inbreeding;
- future: effective founder contribution, genomic relationship matrix.

## L5 — Decision optimization

Chứa:

- target phenotype probability;
- expected information gain;
- uncertainty penalty;
- health/welfare penalties;
- hard relation constraints;
- pair ranking.

L5 không được tự tạo consent hoặc bỏ block.

## L6 — Phenotype/pigment screen

Chứa:

- feature-level QA;
- T02 rule engine;
- artifact/UNK guard;
- future: calibrated image pipeline.

L6 không import genotype inference để tránh circular reasoning.

## L7 — Application service

Chứa:

- request routing;
- orchestration;
- envelope;
- error mapping;
- report generation.

L7 không chứa công thức di truyền mới.

## L8 — Adapters

Hiện tại:

- JSON file;
- CLI.

Tương lai:

- PostgreSQL/Supabase;
- Google Sheet import;
- object storage;
- FastAPI/REST;
- batch worker;
- notebook adapter.

## L9 — Observability và governance

- input/output hash;
- run ID;
- model version;
- configuration snapshot;
- reviewer status;
- incident log;
- data lineage.

---

## 5. Dependency rule

```text
adapters ───────► application ───────► domain/core
                      │                 ▲
                      ├────► inference ─┤
                      ├────► pedigree ──┤
                      ├────► optimizer ─┤
                      └────► pigment ───┘
```

Không cho phép:

- core import CLI;
- inference import database;
- pigment import genotype posterior;
- optimizer sửa pedigree;
- report writer thay numerical result;
- model coder thêm network call trong deterministic function.

---

## 6. Bounded contexts

### 6.1. Phenotype Context

Nguồn sự thật:

- T02 ontology;
- ảnh gốc/Q-level;
- timeline thay lông;
- reviewer agreement.

Output:

- phenotype label probability;
- không có genotype.

### 6.2. Breeding Record Context

Nguồn sự thật:

- pairing ID;
- clutch ID;
- parent ID;
- counts theo mốc;
- sex method;
- survival.

Output:

- observation table cho inference.

### 6.3. Genetics Context

Nguồn sự thật:

- model registry;
- genotype posterior;
- offspring observations;
- molecular result.

Output:

- model support;
- cross distribution;
- evidence grade.

### 6.4. Pedigree Context

Nguồn sự thật:

- sire/dam links;
- founder IDs;
- unknown parent policy.

Output:

- relationship, kinship, inbreeding, relation blocks.

### 6.5. Welfare Context

Nguồn sự thật:

- veterinary review;
- health/fertility flags;
- consent/source status.

Output:

- eligible/blocked.

Không để genetics score override welfare block.

---

## 7. State machines

### 7.1. Case state

```text
received
  ↓
identity_validated
  ↓
phenotype_reviewed
  ↓
genetics_data_pending ───────────────┐
  ↓                                  │
genetics_inference_ready             │
  ↓                                  │
human_review_required                │
  ├──► approved_for_research_cross    │
  ├──► needs_more_data ──────────────┘
  └──► blocked
```

### 7.2. Model lifecycle

```text
draft
  ↓
unit_tested
  ↓
scientifically_reviewed
  ↓
validated_on_holdout
  ↓
active
  ↓
deprecated
```

Mỗi run phải lưu `model_id` + `model_version`.

### 7.3. Evidence lifecycle

```text
G0 → G1 → G2 → G3 → G4 → G5 → G6
```

Không tự động nâng chỉ vì posterior cao. Evidence phụ thuộc loại dữ liệu, không chỉ numerical confidence.

---

## 8. Model registry

Mỗi mô hình cần metadata:

```yaml
model_id: M4_ZR
version: 0.1.0
mode: z_linked_recessive
loci:
  - locus_id: symbolic_M
    chromosome: Z
    alleles: [N, m]
phenotype_map_version: 0.1.0
penetrance: 1.0
observation_error: 0.01
assumptions:
  - one_locus
  - two_alleles
  - equal_transmission
  - no_genotype_specific_survival
validated_species_scope: generic_avian_ZW_demo
status: research_only
```

### 8.1. Registry invariant

- `model_id` định danh family;
- `version` định danh semantics;
- config thay penetrance có thể là run config, nhưng thay phenotype mapping phải tăng model version;
- không overwrite model active; tạo bản mới;
- report phải giữ snapshot model.

---

## 9. Plugin interface mở rộng

Để thêm mô hình mới, production code nên dùng protocol:

```python
class InheritancePlugin(Protocol):
    model_id: str
    version: str

    def enumerate_parent_genotypes(self, sex: Sex) -> list[Genotype]: ...
    def gametes(self, genotype: Genotype, sex: Sex) -> Distribution: ...
    def offspring(self, paternal_gamete, maternal_gamete) -> OffspringState: ...
    def phenotype_distribution(self, genotype, sex, covariates) -> Distribution: ...
    def validate(self) -> list[Issue]: ...
```

Các plugin dự kiến:

- one-locus autosomal;
- one-locus Z-linked;
- two-locus epistasis;
- mitochondrial/maternal effect nếu có dữ liệu;
- polygenic threshold;
- age-dependent penetrance;
- somatic mosaic/non-heritable.

Plugin không được tự gọi optimizer hoặc viết database.

---

## 10. Luồng dữ liệu chi tiết

### 10.1. Cross simulation

```text
request JSON
  ↓ schema/domain validation
model snapshot
  ↓
parent genotype posterior
  ↓ normalize
all genotype-pair combinations
  ↓
gamete distributions
  ↓
offspring genotype distribution
  ↓
phenotype distribution + penetrance
  ↓
sex-specific output
  ↓ evidence/limits envelope
```

### 10.2. Inference

```text
parent visible phenotype
offspring counts
candidate models + priors
          │
          ▼
compatible parent genotypes per model
          │
          ▼
all parent configurations
          │
          ▼
predicted offspring categories
          │
          ▼
observation-error smoothing
          │
          ▼
multinomial likelihood
          │
          ▼
marginalize parent configuration
          │
          ▼
model posterior/ranking
```

### 10.3. Mate ranking

```text
model weights + genotype posteriors
              │
              ▼
   predicted distribution per model
              │
     ┌────────┼─────────┐
     ▼        ▼         ▼
 target P    EIG     uncertainty
     │        │         │
     └────────┼─────────┘
              ▼
 pedigree + health/welfare constraints
              ▼
 multi-objective score
              ▼
 eligible ranking + blocked list
```

---

## 11. Dữ liệu và persistence production

### 11.1. Bảng đề xuất

```text
cases
individuals
phenotype_observations
media_assets
pedigree_links
pairings
clutches
offspring_observations
inheritance_models
model_runs
genotype_posteriors
mate_ranking_runs
molecular_samples
molecular_results
reviews
incidents
```

### 11.2. Event sourcing tối thiểu

Không sửa đè các sự kiện:

- phenotype reclassification;
- sex correction;
- pedigree correction;
- offspring count correction;
- model version update;
- consent withdrawal.

Lưu event:

```yaml
event_id:
entity_type:
entity_id:
event_type:
previous_version:
new_version:
reason:
actor:
timestamp:
input_hash:
```

### 11.3. Immutable run

Một `model_run` phải giữ:

- input snapshot/hash;
- code commit/version;
- model config;
- seed nếu có Monte Carlo;
- output;
- warnings;
- reviewer status.

---

## 12. API production đề xuất

### 12.1. Endpoints

```text
POST /v1/cross/simulate
POST /v1/inheritance/infer
POST /v1/pedigree/evaluate
POST /v1/mates/rank
POST /v1/phenotype/rank
GET  /v1/models
GET  /v1/runs/{run_id}
POST /v1/runs/{run_id}/review
```

### 12.2. Idempotency

- POST nhận `idempotency_key`;
- cùng key + cùng input hash trả run cũ;
- cùng key + input khác trả conflict;
- không tạo hai model run cho cùng immutable request.

### 12.3. Status codes logic

- `200`: deterministic result;
- `202`: batch/long-running accepted;
- `400`: schema invalid;
- `409`: version/idempotency conflict;
- `422`: biologically inconsistent domain input;
- `423`: blocked by consent/welfare;
- `500`: unexpected internal error, không trả partial claim.

---

## 13. Error taxonomy

```text
input_json_invalid
schema_version_unsupported
individual_id_missing
duplicate_offspring_id
sex_required_for_z_linked_model
parent_phenotype_incompatible_with_model
pedigree_cycle_detected
pedigree_individual_unknown
full_sibling_pair
direct_ancestor_pair
offspring_inbreeding_above_configured_limit
missing_genotype_posteriors
observation_categories_empty
model_priors_zero_total
image_quality_too_low
claim_not_allowed
```

Error phải machine-readable và có human message riêng ở UI; core không hard-code lời khuyên thương mại.

---

## 14. Bất định và calibration

### 14.1. Ba loại bất định

- `aleatoric`: ngẫu nhiên phân ly, survival, phenotype noise;
- `epistemic`: chưa biết genotype/model;
- `data_quality`: thiếu giới, ảnh Q thấp, pedigree thiếu.

Không gộp ba loại thành một confidence duy nhất.

### 14.2. Output đề xuất

```yaml
uncertainty:
  segregation_randomness: described
  parent_genotype_entropy: 0.42
  model_entropy: 1.12
  phenotype_measurement_quality: Q2
  missingness_flags: [sex_unknown, first_molt_missing]
```

### 14.3. Calibration

Khi có đủ case:

- reliability diagram cho posterior model;
- Brier score cho phenotype prediction;
- coverage–accuracy cho abstention;
- external validation theo cơ sở nuôi/khu vực/thiết bị;
- không chia ảnh cùng cá thể vào train và test.

---

## 15. Security, privacy và misuse

- Tách PII người nuôi khỏi ID chim.
- Không công khai địa chỉ chính xác cá thể hiếm.
- Không đưa ảnh vào training nếu media licence không cho phép.
- Không nhận source claim là legal verification.
- Không để user sửa `welfare_blocked:false` nếu cờ đến từ reviewer có thẩm quyền; production cần role-based override.
- Log mọi lần bỏ block và lý do.
- Không để LLM tạo SQL trực tiếp vào bảng production nếu chưa qua allowlisted service.
- Không dùng service credential trong skill docs hoặc output.

---

## 16. Testing architecture

### 16.1. Unit tests

- genotype parser;
- gamete distribution;
- Punnett expected ratios;
- Z-linked reciprocal crosses;
- penetrance;
- likelihood;
- pedigree matrix;
- EIG;
- evidence guard.

### 16.2. Property tests tương lai

- tổng xác suất = 1;
- đổi thứ tự allele không đổi kết quả;
- cross autosomal đổi cha/mẹ không đổi distribution nếu không có sex effect;
- Z-linked reciprocal cross có thể khác;
- posterior model sum = 1;
- F đời con nằm `[0,1]`;
- hard-block luôn thắng score.

### 16.3. Golden tests

Lưu các case chuẩn:

- AR `N/m × m/m`;
- ZR `m/m × N/W`;
- reciprocal ZR;
- sex-skewed brood;
- full-sibling block;
- pigment PGR;
- low-Q image → UNK.

### 16.4. Mutation testing tương lai

Đặc biệt kiểm:

- đảo male/female trong Z-linked;
- bỏ hệ số 0.5 trong relationship matrix;
- bỏ observation error;
- cho block pair vẫn đứng hạng eligible.

---

## 17. Performance

Prototype:

- một locus: O(Gm × Gf × gametes);
- inference: O(models × compatible parent configurations × categories);
- pedigree matrix: O(n²);
- pair ranking: O(males × females × models × genotype combinations).

Khi lớn:

- cache cross distribution theo `(model_version, male_genotype, female_genotype)`;
- sparse pedigree hoặc specialized library;
- batch vectorization;
- queue cho rank hàng nghìn cặp;
- precompute model predictions;
- không cache kết quả nếu consent/version đã thay.

---

## 18. Evolution path

### V0.1 — Prototype deterministic

- hoàn thành trong package này;
- local JSON;
- standard library;
- one-locus.

### V0.2 — Data validation và two-locus

- robust schema validator;
- epistasis plugin;
- survival/missingness fields;
- richer evidence report.

### V0.3 — Database/API

- Postgres schema;
- FastAPI adapter;
- run registry;
- role/consent gates;
- immutable audit.

### V0.4 — Image feature service

- RAW/calibration ingestion;
- bird/body segmentation;
- CIELAB/spectral feature extraction;
- out-of-distribution guard;
- reviewer labeling UI.

### V0.5 — Molecular integration

- sample/chain-of-custody;
- variant table;
- cosegregation;
- genomic relationship matrix;
- G5 report.

### V1.0 — Validated research platform

Chỉ đạt khi:

- external validation;
- model registry/governance;
- reproducible deployment;
- privacy/consent audit;
- veterinary/scientific review;
- documented failure modes.

---

## 19. Quyết định thiết kế cốt lõi

1. **Core dùng ký hiệu allele, không tên gene.** Giảm nguy cơ biến giả thuyết thành fact.
2. **Posterior genotype là first-class object.** Không ép carrier khi chưa có DNA/test cross.
3. **M0_NONGENETIC luôn tồn tại.** Tránh gene-first bias.
4. **Phenotype engine tách khỏi genetics engine.** Tránh circular evidence.
5. **Hard constraints trước score.** Cặp bị block không thể thắng nhờ màu.
6. **Evidence grade độc lập với numerical confidence.** Posterior cao không tạo G5.
7. **One-file prototype + modular package.** Một bản dễ đọc, một bản dễ mở rộng.
8. **Không dependency ngoài ở v0.1.0.** Dễ chạy, audit và chuyển giao.

---

## 20. Definition of architecture done

Một implementation được coi là phù hợp kiến trúc khi:

- giữ dependency rule;
- có schema/version;
- deterministic core không network;
- model output có limits/evidence;
- Z-linked bắt buộc giới;
- pedigree cycle bị chặn;
- direct/full-sibling block không bị score override;
- image/phenotype không tạo genotype;
- test toàn bộ invariant;
- thay semantics có migration/changelog.


<!-- END ARCHITECTURE -->

---

<!-- BEGIN SCIENTIFIC BOUNDARIES: references/scientific-boundaries.md -->

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

Skill v0.1.0 luôn trả `biochemical_pigment_status:unconfirmed`.

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


<!-- END SCIENTIFIC BOUNDARIES -->

---

<!-- BEGIN DATA CONTRACTS: references/data-contracts.md -->

# Data contracts — BioBirdDNA Genetics

Mục tiêu của tài liệu này là để model/coder khác tạo đúng input/output mà không phải đoán semantics.

---

## 1. Versioning

Mỗi request production nên có:

```json
{
  "schema_version": "1.0.0",
  "request_id": "CBIM-GEN-RUN-...",
  "case_id": "CM-CASE-..."
}
```

Prototype v0.1.0 chưa bắt buộc hai field này trong CLI để giữ ví dụ gọn, nhưng production phải thêm.

### SemVer

- `PATCH`: bug fix không đổi semantics.
- `MINOR`: thêm field optional/mô hình mới.
- `MAJOR`: đổi nghĩa field hoặc output breaking.

---

## 2. ID conventions

```text
Case:       CM-CASE-000001
Individual: CM-IND-000001
Pairing:    CM-MAT-000001
Clutch:     CM-CLU-000001
Offspring:  CM-OFF-000001
Sample:     CM-SMP-000001
Model run:  CM-RUN-000001
Review:     CM-REV-000001
```

ID không chứa tên chủ chim, số điện thoại hoặc địa chỉ.

---

## 3. Inheritance model

```json
{
  "model_id": "M4_ZR",
  "mode": "z_linked_recessive",
  "normal_allele": "N",
  "variant_allele": "m",
  "normal_label": "normal",
  "target_label": "target_colour",
  "intermediate_label": "intermediate",
  "penetrance": 1.0,
  "observation_error": 0.01,
  "female_probability": 0.5,
  "baseline_target_probability": 0.05,
  "description": "..."
}
```

### `mode`

```text
autosomal_recessive
autosomal_dominant
incomplete_dominant
z_linked_recessive
z_linked_dominant
non_genetic
```

### Invariant

- `normal_allele != variant_allele`.
- probability fields trong `[0,1]`.
- `model_id` + model version production phải immutable.

---

## 4. Genotype strings

### Autosomal hoặc trống ZZ

```text
N/N
N/m
m/m
```

Parser canonicalize `m/N` thành `N/m`.

### Mái ZW

```text
N/W
m/W
```

Không dùng `W/N` hoặc `W/m`.

### Non-genetic

```text
NA
```

---

## 5. Cross request

### Known genotypes

```json
{
  "model": {"model_id": "M1_AR", "mode": "autosomal_recessive"},
  "male_genotype": "N/m",
  "female_genotype": "m/m"
}
```

### Uncertain genotypes

```json
{
  "model": {"model_id": "M1_AR", "mode": "autosomal_recessive"},
  "male_genotypes": {"N/N": 0.3, "N/m": 0.7},
  "female_genotypes": {"N/m": 1.0}
}
```

Distribution không cần sum chính xác 1; engine normalize. Tổng phải >0 và không có số âm.

---

## 6. Cross output

```json
{
  "model_id": "M1_AR",
  "mode": "autosomal_recessive",
  "male_genotype": "N/m",
  "female_genotype": "m/m",
  "genotype_probabilities": {
    "male:N/m": 0.25,
    "male:m/m": 0.25,
    "female:N/m": 0.25,
    "female:m/m": 0.25
  },
  "phenotype_probabilities": {
    "male:normal": 0.25,
    "male:target_colour": 0.25,
    "female:normal": 0.25,
    "female:target_colour": 0.25
  },
  "limits": {
    "gene_status": "model_assumption_only",
    "inheritance_status": "simulation_not_confirmation"
  }
}
```

Keys có dạng:

```text
<sex>:<genotype>
<sex>:<phenotype>
```

---

## 7. Inference request

```json
{
  "target_label": "target_colour",
  "male_parent_phenotype": "target_colour",
  "female_parent_phenotype": "normal",
  "offspring": [
    {"sex": "male", "phenotype": "normal", "count": 5},
    {"sex": "male", "phenotype": "target_colour", "count": 0},
    {"sex": "female", "phenotype": "normal", "count": 0},
    {"sex": "female", "phenotype": "target_colour", "count": 6}
  ]
}
```

### Sex mode

Một request phải dùng một trong hai kiểu:

- tất cả `male/female`;
- tất cả `unknown`.

Không trộn known và unknown sex trong một request v0.1.0.

### Counts

- integer ≥0;
- zero là quan sát thật rằng category đã được theo dõi nhưng không xuất hiện;
- field bị thiếu nghĩa là category chưa được khai báo, không tự động giống zero trong mọi analysis production.

---

## 8. Inference output

```json
{
  "ranking": [
    {
      "model_id": "M4_ZR",
      "mode": "z_linked_recessive",
      "prior": 0.1667,
      "posterior": 0.98,
      "log_evidence": -7.1,
      "best_parent_configuration": {
        "male_genotype": "m/m",
        "female_genotype": "N/W",
        "posterior_within_model": 0.99
      }
    }
  ],
  "best_model": "M4_ZR",
  "best_model_posterior": 0.98,
  "limits": {
    "claim": "candidate_inheritance_model_only",
    "gene_status": "not_assessed"
  }
}
```

`posterior` không được đổi tên thành `gene_probability`.

---

## 9. Pedigree object

```json
{
  "individual_id": "CM-IND-000001",
  "sex": "male",
  "sire_id": "CM-IND-000010",
  "dam_id": "CM-IND-000011"
}
```

Unknown parent dùng `null` hoặc bỏ field. Production nên có thêm:

```text
parentage_method
parentage_confidence
parentage_source
```

---

## 10. Genetic profile

```json
{
  "individual_id": "CM-IND-000001",
  "sex": "male",
  "phenotype_label": "normal",
  "health_risk": 0.05,
  "fertility_risk": 0.0,
  "welfare_blocked": false,
  "genotype_posteriors": {
    "M1_AR": {"N/N": 0.2, "N/m": 0.8},
    "M4_ZR": {"N/N": 0.3, "N/m": 0.7}
  }
}
```

### `health_risk`/`fertility_risk`

Prototype nhận số `[0,1]`, nhưng production phải lưu:

- source;
- reviewer;
- date;
- rationale;
- expiry.

Không cho LLM tự tạo risk từ mô tả mơ hồ.

---

## 11. Mate ranking request

Các khối:

```text
pedigree[]
males[]
females[]
hypotheses[]
weights{}
max_offspring_inbreeding
block_direct_ancestor_and_full_sibling
```

Mỗi hypothesis:

```json
{
  "weight": 0.6,
  "model": {"model_id": "M1_AR", "mode": "autosomal_recessive"}
}
```

Weights model sẽ được normalize.

---

## 12. Mate ranking output

```json
{
  "pair_ranking": [
    {
      "male_id": "M-001",
      "female_id": "F-002",
      "eligible": true,
      "block_reasons": [],
      "score": 0.31,
      "target_probability": 0.28,
      "expected_information_gain_bits": 0.14,
      "offspring_inbreeding": 0.0,
      "health_risk": 0.05,
      "genotype_uncertainty": 0.63
    }
  ],
  "limits": {
    "purpose": "decision_support_not_automatic_breeding_order"
  }
}
```

Blocked pair vẫn xuất trong list.

---

## 13. Pigment request

Tất cả feature `[0,1]` trừ bool/null/Q-level.

```json
{
  "white_area_ratio": 0.18,
  "brown_shift_score": 0.02,
  "grey_shift_score": 0.03,
  "dark_expansion_score": 0.0,
  "red_patch_change_score": 0.0,
  "pattern_retained": false,
  "iris_red_confirmed": false,
  "flash_red_eye_risk": false,
  "progressive_across_molts": true,
  "present_in_first_plumage": false,
  "localized_after_injury": false,
  "bilateral_symmetry": 0.35,
  "image_quality": "Q3"
}
```

Các score này phải đến từ pipeline đo/reviewer, không tự suy từ caption.

---

## 14. Pigment output

```json
{
  "primary_hypothesis": "MEL-PGR",
  "alternative_hypotheses": ["MEL-LEU", "UNK"],
  "scores": {"MEL-PGR": 0.55, "MEL-LEU": 0.22},
  "rationale": {"MEL-PGR": ["white_or_pale_area_progresses_across_molts"]},
  "limits": {
    "measurement_domain": "visible_phenotype_only",
    "biochemical_pigment_status": "unconfirmed",
    "gene_status": "not_assessed",
    "inheritance_status": "not_assessed"
  }
}
```

`score` là normalized heuristic, không phải gene posterior.

---

## 15. Evidence object production

```json
{
  "phenotype_evidence": "E2",
  "genetics_evidence": "G2",
  "image_quality": "Q3",
  "confidence": "C2",
  "claim_ceiling": "candidate_inheritance_model",
  "rationale": ["two_repeated_broods_same_pair"],
  "missing_for_next_level": ["independent_pair", "reciprocal_cross"]
}
```

---

## 16. Missing data contract

Không dùng empty string cho unknown. Dùng:

- `null` khi unknown;
- enum `unknown` khi cần giữ explicit category;
- `missing_data[]` giải thích field nào và ảnh hưởng gì.

Ví dụ:

```json
{
  "missing_data": [
    {
      "field": "offspring_sex",
      "impact": "cannot_discriminate_Z_linkage",
      "priority": "high"
    }
  ]
}
```

---

## 17. Run envelope production

```json
{
  "ok": true,
  "operation": "infer",
  "engine_version": "0.1.0",
  "schema_version": "1.0.0",
  "run_id": "CM-RUN-...",
  "input_sha256": "...",
  "model_snapshots": [],
  "result": {},
  "warnings": [],
  "human_review_status": "pending"
}
```

---

## 18. JSON Schemas

Xem:

- `assets/schemas/inheritance-model.schema.json`
- `assets/schemas/cross-request.schema.json`
- `assets/schemas/inference-request.schema.json`
- `assets/schemas/mate-ranking-request.schema.json`
- `assets/schemas/pigment-request.schema.json`

Prototype không phụ thuộc thư viện JSON Schema; schema hiện phục vụ contract, validation ngoài và model handoff.


<!-- END DATA CONTRACTS -->

---

<!-- BEGIN ALGORITHMS: references/algorithms.md -->

# Thuật toán và công thức — BioBirdDNA Genetics

Tài liệu này mô tả lõi toán học để model khác có thể tái triển khai mà không làm thay đổi semantics.

---

## 1. Quy ước

- Chim trống: `ZZ`.
- Chim mái: `ZW`.
- `N`: allele ký hiệu bình thường.
- `m`: allele ký hiệu biến đổi.
- `M`: mô hình di truyền.
- `G♂`, `G♀`: genotype cha/mẹ.
- `Y`: phenotype đời con.
- `D`: dữ liệu quan sát.

`N` và `m` không phải tên gene thật.

---

## 2. Phân bố giao tử

Với genotype autosomal hoặc ZZ:

```text
N/N → P(N)=1
N/m → P(N)=0.5, P(m)=0.5
m/m → P(m)=1
```

Với mái ZW:

```text
N/W → P(ZN)=0.5, P(W)=0.5
m/W → P(Zm)=0.5, P(W)=0.5
```

Prototype giả định segregation đều. Nếu có meiotic drive hoặc transmission distortion, phải thêm plugin khác và không sửa ngầm engine hiện tại.

---

## 3. Cross genotype chắc chắn

### 3.1. NST thường


a) Chọn giao tử từ cha và mẹ.  
 b) Kết hợp thành genotype đời con.  
 c) Nhân thêm xác suất giới.

\[
P(G_c, S_c) = P(a_{cha}) P(a_{mẹ}) P(S_c)
\]

Mặc định:

\[
P(S_c=mái)=P(S_c=trống)=0.5
\]

### 3.2. Liên kết Z

- Mẹ truyền `W` → con mái; Z của con mái đến từ cha.
- Mẹ truyền `Z` → con trống; con trống nhận một Z từ mỗi parent.

\[
P(G_c,S_c)=P(Z_{cha})P(Z/W_{mẹ})
\]

Không nhân sex ratio riêng vì sex đã phát sinh từ giao tử mẹ.

---

## 4. Genotype → phenotype

### 4.1. Lặn thường

```text
N/N → normal
N/m → normal visible; carrier state
m/m → target
```

### 4.2. Trội thường

```text
N/N → normal
N/m → target
m/m → target
```

### 4.3. Trội không hoàn toàn

```text
N/N → normal
N/m → intermediate
m/m → target
```

### 4.4. Lặn liên kết Z

Trống:

```text
N/N → normal
N/m → normal visible; carrier state
m/m → target
```

Mái:

```text
N/W → normal
m/W → target
```

### 4.5. Trội liên kết Z

Trống:

```text
N/N → normal
N/m hoặc m/m → target
```

Mái:

```text
N/W → normal
m/W → target
```

---

## 5. Penetrance

Nếu genotype thuộc state affected nhưng penetrance là `π`:

\[
P(Y=target|G)=\pi
\]

\[
P(Y=normal|G)=1-\pi
\]

Prototype dùng một penetrance chung. Production nên cho phép penetrance theo:

- genotype;
- giới;
- tuổi;
- mùa lông;
- môi trường;
- dòng gia đình.

Không dùng penetrance để che lấp chấm phenotype kém; hai nguồn bất định phải tách.

---

## 6. Parent genotype chưa chắc chắn

Khi cha/mẹ có posterior:

\[
P(G_{cha}|D_{cha}), P(G_{mẹ}|D_{mẹ})
\]

phân bố đời con:

\[
P(Y)=\sum_{G_{cha},G_{mẹ},G_c}
P(Y|G_c,M)
P(G_c|G_{cha},G_{mẹ},M)
P(G_{cha}|D_{cha})
P(G_{mẹ}|D_{mẹ})
\]

Code phải normalize posterior đầu vào và lỗi nếu tổng bằng 0.

---

## 7. Likelihood đời con

Với các category `k=1..K`, số quan sát `n_k`, xác suất model dự đoán `p_k`:

\[
P(D|M,G_{cha},G_{mẹ}) =
\frac{N!}{\prod_k n_k!}\prod_k p_k^{n_k}
\]

Log-likelihood:

\[
\log L = \log N! - \sum_k \log n_k! + \sum_k n_k \log p_k
\]

Code dùng `lgamma(n+1)` để tránh overflow.

---

## 8. Observation error

Để tránh một lỗi chấm duy nhất làm probability bằng 0 tuyệt đối, prototype dùng:

\[
\tilde{p}_k = (1-\epsilon)p_k + \frac{\epsilon}{K}
\]

Sau đó normalize lại trên category set.

`ε` không được đặt lớn để cứu model sai. Nó phải được version/config và sau này estimate từ reviewer disagreement.

---

## 9. Marginalize parent configuration

Nếu có nhiều genotype cha/mẹ tương thích phenotype:

\[
P(D|M)=\sum_c P(D|M,c)P(c|M)
\]

Trong log-space:

\[
\log P(D|M)=\operatorname{logsumexp}
(\log P(D|M,c)+\log P(c|M))
\]

Prototype dùng prior đều giữa configuration khi không có prior riêng.

---

## 10. Model posterior

\[
P(M|D)=\frac{P(D|M)P(M)}{\sum_jP(D|M_j)P(M_j)}
\]

Report phải nêu rõ:

- model set được so sánh;
- prior;
- observation error;
- parent phenotype assumptions;
- dữ liệu bị thiếu.

Posterior chỉ tương đối trong model set. Nếu model thật không nằm trong tập, posterior cao vẫn có thể gây hiểu nhầm.

---

## 11. Non-genetic baseline

Prototype `M0_NONGENETIC` dùng target probability nền `p0` độc lập parent/genotype:

\[
P(Y=target)=p_0
\]

Đây chưa phải model sinh học hoàn chỉnh. Nó là guard để model comparison không bắt buộc mọi data phải mang tính Mendel.

Production nên thay bằng các candidate:

- age-dependent;
- injury-localized;
- progressive greying;
- environment-associated;
- random classification error.

---

## 12. Relationship matrix

Numerator relationship matrix `A`:

Với founder không biết parent:

\[
A_{ii}=1
\]

Với cá thể `i` có sire `s`, dam `d`:

\[
A_{ij}=\frac{1}{2}(A_{sj}+A_{dj}), \quad j<i
\]

\[
A_{ii}=1+\frac{1}{2}A_{sd}
\]

Nếu một parent không biết, contribution parent đó là 0 trong prototype.

### 12.1. Inbreeding cá thể

\[
F_i=A_{ii}-1
\]

### 12.2. Inbreeding đời con dự kiến

\[
F_{con}=\frac{1}{2}A_{cha,mẹ}
\]

Founders “không biết” được coi không liên quan là giả định tiện dụng, có thể đánh giá thấp cận huyết. Report production phải gắn flag `founder_relatedness_unknown`.

---

## 13. Expected information gain

Với model hypothesis `M` có weight `w_m` và outcome distribution `P(Y|M,pair)`:

Mixture:

\[
P(Y)=\sum_m w_mP(Y|M,pair)
\]

Expected information gain cho một offspring:

\[
EIG = H(P(Y)) - \sum_m w_mH(P(Y|M,pair))
\]

Trong đó:

\[
H(P)=-\sum_yP(y)\log_2P(y)
\]

EIG cao khi các model dự đoán outcome khác nhau. Nó không đồng nghĩa cặp đó tốt về sức khỏe hoặc đạo đức.

---

## 14. Genotype uncertainty

Với posterior genotype `q_g`:

\[
H_G=-\sum_g q_g\log_2q_g
\]

Normalize:

\[
U_G=\frac{H_G}{\log_2 |G|}
\]

`U_G` nằm `[0,1]` khi có trên một genotype state.

---

## 15. Pair score

Prototype:

\[
Score =
 w_TP(target)
 + w_IEIG
 - w_FF_{con}
 - w_HRisk
 - w_UU_G
\]

Thứ tự xử lý:

1. hard block relation/welfare;
2. tính score;
3. sort eligible trước;
4. giữ blocked pair trong report.

Không được “phạt nhẹ” cha–con thay vì block.

---

## 16. Rule engine sắc tố

Feature đầu vào là số đã đo hoặc reviewer coding:

- white area ratio;
- brown/grey/dark shift;
- red-domain change;
- pattern retained;
- iris red confirmed;
- flash risk;
- progressive across molts;
- first plumage;
- injury localization;
- Q-level.

Rule cộng score rồi normalize:

\[
P_{heuristic}(label)=\frac{s_{label}}{\sum_js_j}
\]

Đây không phải posterior sinh học. Tên trường `scores` được dùng thay vì `probability_of_gene`.

---

## 17. Monte Carlo tương lai

Analytic cross đủ cho one-locus. Monte Carlo chỉ cần khi:

- multi-locus lớn;
- uncertain survival;
- family-level random effects;
- complex age penetrance;
- distribution of brood outcomes.

Mỗi run phải lưu seed. Không dùng Monte Carlo nếu analytic exact có thể tính.

---

## 18. Confidence interval cho một lứa tương lai

Nếu target probability là `p`, số offspring `n`:

\[
K\sim Binomial(n,p)
\]

Có thể xuất distribution `P(K=k)`, không hứa “n × p con”.

Khi nhiều phenotype category:

\[
\mathbf{K}\sim Multinomial(n,\mathbf{p})
\]

---

## 19. Survival model tương lai

Cần tách:

\[
P(phenotype\ observed)=
P(genotype\ at\ conception)
\times P(survival|genotype)
\times P(classified|phenotype)
\]

Nếu chỉ dùng chim sống sau thay lông, không thể mặc định tỷ lệ đó bằng segregation lúc thụ tinh.

---

## 20. Acceptance tests toán học

- mọi distribution sum gần 1;
- AR `N/m × m/m` → target 0.5;
- ZR `m/m × N/W` → female target 0.5, male normal 0.5;
- full siblings A=0.5, offspring F=0.25;
- unrelated founders F=0;
- posterior models sum 1;
- evidence guard chặn gene claim ở G4.


<!-- END ALGORITHMS -->

---

<!-- BEGIN T02 INTEGRATION: references/t02-integration.md -->

# Tích hợp T02 vào BioBirdDNA Genetics Skill

Nguồn máy đọc được bundle tại:

```text
knowledge/T02_knowledge_base_v1.0.json
```

T02 định nghĩa 14 lớp ontology, thang E/C/Q, trường phenotype, timeline và nguyên tắc không suy genotype từ ảnh.

---

## 1. Mapping lớp T02

Các nhãn phenotype screen v0.1.0:

```text
MEL-LEU
MEL-PGR
MEL-ALB
MEL-INO
MEL-BRN
MEL-DIL
MEL-MEL
CAR-ALT
LOC-ACQ
IMG-ART
BASE-NORM
UNK
```

`PAT-MOS` và `TAX-UNK` nằm trong T02 nhưng chưa có rule feature riêng trong prototype; production phải thêm trước khi tuyên bố coverage toàn ontology.

---

## 2. Trường T02 cần cho genetics

Tối thiểu:

```text
individual_id
sex
sex_method
age_class
molt_stage
primary_hypothesis
alternative_hypotheses
evidence_level
confidence_level
image_quality_score
sire_id
dam_id
mating_records
offspring_counts
```

Nhưng inference production nên lấy dữ liệu từ bảng breeding records chuẩn, không parse text `mating_records`.

---

## 3. Handoff contract

T02 → Genetics:

```json
{
  "individual_id": "CM-IND-...",
  "phenotype_label": "MEL-PGR",
  "phenotype_label_probabilities": {
    "MEL-PGR": 0.55,
    "MEL-LEU": 0.25,
    "UNK": 0.20
  },
  "image_quality": "Q3",
  "phenotype_evidence": "E2",
  "timeline": {
    "progressive_across_molts": true
  }
}
```

Genetics không được đổi `MEL-PGR` thành allele `m` tự động. Cần một model hypothesis riêng nói phenotype target nào đang được theo dõi.

---

## 4. Community name

Giữ:

```text
community_name_verbatim
```

nhưng không dùng làm feature genotype inference. Nó chỉ dùng:

- tìm kiếm;
- hiển thị;
- nghiên cứu ngôn ngữ;
- alias mapping.

---

## 5. Quality ceiling

- Q0–Q1: không measurement-ready; ưu tiên `UNK`.
- Q2: screening thị giác.
- Q3: định lượng ảnh chuẩn hơn.
- Q4: calibrated research image.

Q-level không thay G-level.

---

## 6. Evidence ceiling

T02 E-level mô tả phenotype/data. Genetics G-level mô tả bằng chứng di truyền. Giữ hai scale riêng:

```text
E2 + G0 có thể xảy ra: phenotype theo dõi tốt nhưng chưa có breeding data.
E1 + G3 có thể bất thường và cần review: breeding data mạnh nhưng phenotype coding ảnh yếu.
```

Final claim ceiling lấy mắt xích yếu nhất liên quan claim.

---

## 7. Model output contract

T02 yêu cầu:

- observations trước assessment;
- primary + alternatives;
- evidence/confidence/quality;
- missing data;
- safe summary.

Genetics report phải giữ cấu trúc tương tự và thêm:

- parent genotype posterior;
- candidate inheritance models;
- offspring distribution;
- pedigree risk;
- genetics grade G0–G6.


<!-- END T02 INTEGRATION -->

---

<!-- BEGIN IMPLEMENTATION ROADMAP: references/implementation-roadmap.md -->

# Lộ trình triển khai production

Tài liệu này biến kiến trúc thành các pha mà model/coder khác có thể thực hiện tuần tự.

---

## Pha 0 — Khóa invariant và test nền

### Mục tiêu

Bảo đảm mọi model khác hiểu cùng một semantics trước khi thêm tính năng.

### Việc làm

- Đọc `SKILL.md` và `scientific-boundaries.md`.
- Chạy toàn bộ unit test.
- Tạo golden output cho năm example.
- Ghi SHA-256 manifest.
- Không thêm framework/network.

### DoD

- test sạch;
- demo chạy;
- không có claim vượt evidence;
- output sum probability = 1.

---

## Pha 1 — Harden deterministic core

### Scope

- robust genotype object thay string nội bộ;
- exact rational/Fraction option cho Punnett;
- custom sex ratio autosomal;
- penetrance theo genotype/sex;
- phenotype state tách carrier visible status;
- property tests.

### Không làm

- database;
- image ML;
- multi-locus.

### DoD

- backward compatible JSON;
- benchmark 10.000 cross;
- no floating drift > tolerance;
- mutation tests cho Z-linked.

---

## Pha 2 — Breeding records và inference v0.2

### Scope

- `pairings`, `clutches`, `offspring` schema;
- sex method/confidence;
- multiple broods;
- survival stages;
- missing category semantics;
- configurable model priors;
- Bayes factor/log evidence report;
- parent genotype posterior export.

### DoD

- ingest ít nhất 20 synthetic pedigrees;
- recover known generating model trong simulation study;
- calibration report;
- abstain khi model set misspecified.

---

## Pha 3 — Multi-locus plugin

### Scope

- generic locus object;
- independent assortment;
- two-locus epistasis phenotype map;
- optional linkage/recombination fraction;
- plugin registry;
- JSON model DSL.

### DoD

- reproduce known synthetic ratios;
- reject invalid phenotype maps;
- complexity guard;
- no gene-name claim.

---

## Pha 4 — Pedigree/population governance

### Scope

- unknown-founder uncertainty;
- half-sibling/first-cousin labels;
- founder contribution;
- mean kinship;
- effective population size proxy;
- limit repeated parent contribution;
- genomic relationship adapter future interface.

### DoD

- test pedigrees with loops rejected;
- relationship matrix compared with independent implementation;
- hard block cannot be overridden by weights.

---

## Pha 5 — Mate optimizer v0.2

### Scope

- separate research and production objectives;
- Pareto frontier instead of only scalar score;
- multi-generation simulation;
- diversity constraints;
- capacity constraints per breeder;
- explainable ranking.

### DoD

- output eligible, blocked, Pareto set;
- EIG validated by simulation;
- report sensitivity to weights;
- no auto-execution.

---

## Pha 6 — Image phenotype pipeline

### Scope

- media provenance and consent gate;
- RAW/JPEG QA;
- calibration target detection;
- bird segmentation;
- body region segmentation;
- color feature extraction;
- longitudinal alignment;
- OOD/UNK;
- two-reviewer labeling.

### Dependency candidates

Chọn sau technical review; không khóa ở blueprint. Core phải nhận feature chuẩn, không phụ thuộc CV library cụ thể.

### DoD

- split theo individual/capture session;
- external device validation;
- feature repeatability;
- model card;
- no gene output.

---

## Pha 7 — API/database

### Scope

- Postgres/Supabase schema;
- migration;
- API adapter;
- auth/RBAC;
- idempotency;
- immutable run registry;
- audit event;
- job queue.

### DoD

- contract tests;
- role tests;
- consent withdrawal path;
- deterministic core remains network-free;
- exact run reproducibility.

---

## Pha 8 — Molecular integration

### Scope

- sample/aliquot IDs;
- chain of custody;
- variant import;
- genotype QC;
- cosegregation;
- family-based association;
- genomic kinship;
- G5 evidence report.

### DoD

- sample swap checks;
- missing genotype handling;
- no G5 without lineage/QC;
- publication language guard.

---

## Pha 9 — Governance và v1.0

### Scope

- scientific review board;
- veterinary review;
- privacy/legal review;
- model registry lifecycle;
- incident response;
- external replication;
- public documentation.

### DoD

- validated dataset;
- model cards;
- data sheets;
- audit trail;
- reproducible deployment;
- documented known failure modes.

---

## Ưu tiên ngay sau v0.1.0

1. Chuẩn hóa breeding record schema.
2. Thêm simulation study cho inference recovery.
3. Tách research pairing và target pairing trong API.
4. Thêm Pareto ranking.
5. Tạo importer Google Sheet/Supabase chỉ sau khi core ổn định.


<!-- END IMPLEMENTATION ROADMAP -->

---

<!-- BEGIN MODEL CARD: references/model-card.md -->

# Model card — BioBirdDNA Genetics v0.1.0

## Intended use

- Educational/research cross simulation.
- Candidate inheritance-model ranking from complete brood data.
- Pedigree and inbreeding decision support.
- Research-pair ranking.
- T02 phenotype screening from precomputed features.

## Out-of-scope use

- Confirming a gene from an image.
- Veterinary diagnosis.
- Automatic breeding order.
- Commercial guarantee of colour outcomes.
- Invasive sample collection guidance.
- Population-wide frequency estimation from social media images.

## Implemented models

- One-locus autosomal recessive.
- One-locus autosomal dominant.
- One-locus incomplete dominance.
- One-locus avian Z-linked recessive.
- One-locus avian Z-linked dominant.
- Simple non-genetic baseline.

## Key assumptions

- Two symbolic alleles.
- Equal transmission.
- No linkage for one-locus models.
- No genotype-specific survival in v0.1.0.
- Simple penetrance.
- Simple observation-error smoothing.
- Unknown founders treated as unrelated in pedigree matrix; this can underestimate true relatedness.

## Known limitations

- Posterior depends on candidate model set and priors.
- Small brood sizes produce unstable ranking.
- Phenotype labels can be age/molt dependent.
- Sex misclassification strongly affects Z-linked inference.
- Non-genetic baseline is intentionally simple.
- Pair score is a configurable decision-support function, not a welfare standard.

## Required human review

- Scientific review for inheritance claims.
- Veterinary/welfare review for real pairings.
- Legal/consent review for data and samples.
- Lab review for molecular claims.

## Evidence language

- G0–G2: exploratory/candidate.
- G3–G4: supported/strong model, no confirmed locus.
- G5: locus/variant supported by cosegregation.
- G6: independently replicated or functionally supported mechanism.


<!-- END MODEL CARD -->

---

<!-- BEGIN MODEL CODER HANDOFF: planning/MODEL_CODER_HANDOFF.md -->

# Handoff cho model/coder tiếp tục phát triển skill

Tài liệu này là prompt vận hành và hợp đồng bàn giao. Model mới phải đọc đúng thứ tự và không bỏ qua invariant.

---

## 1. Read order bắt buộc

1. `SKILL.md`
2. `references/scientific-boundaries.md`
3. `references/architecture.md`
4. `references/data-contracts.md`
5. `references/algorithms.md`
6. `references/t02-integration.md`
7. `planning/backlog.json`
8. code và test liên quan issue.

Không bắt đầu code chỉ từ mô tả issue.

---

## 2. Prompt chuẩn cho coding model

```text
Dùng skill $biobirddna-genetics tại <skill-root>.

Nhiệm vụ: thực hiện đúng issue <ISSUE_ID> trong planning/backlog.json.

Trước khi sửa:
1. Đọc SKILL.md, scientific-boundaries.md, architecture.md, data-contracts.md và algorithms.md.
2. Chạy toàn bộ test hiện tại và ghi baseline.
3. Liệt kê invariant liên quan và file dự kiến sửa.
4. Không đổi public contract ngoài phạm vi issue.

Trong khi sửa:
- Core deterministic không gọi network/database.
- Không suy gene/genotype từ ảnh.
- Không để score override hard block.
- Không xóa test để pass.
- Thêm test fail trước, sau đó code pass.
- Mọi output mới phải có limits/evidence.

Sau khi sửa:
1. Chạy unit test toàn bộ.
2. Chạy demo liên quan.
3. Chạy scripts/validate_package.py.
4. Cập nhật schema/docs/changelog nếu contract hoặc semantics đổi.
5. Báo file đã sửa, test result, limitations và backlog tiếp theo.
```

---

## 3. Quy tắc patch

### Được phép

- thêm module mới theo dependency rule;
- thêm field optional có default;
- thêm model plugin;
- thêm test/golden case;
- cải thiện numerical stability;
- sửa bug có regression test.

### Cần ADR

- đổi representation genotype;
- thêm framework/dependency ngoài;
- đổi database;
- đổi evidence thresholds;
- đổi formula pair score;
- đổi default inbreeding policy;
- cho phép model tự động quyết định ghép.

### Bị cấm

- xóa `M0_NONGENETIC` để tăng posterior model gene;
- đổi `posterior` thành “độ chắc chắn gene”;
- gộp T02 phenotype score với genotype posterior;
- tự suy sex;
- bỏ full-sibling/direct-ancestor block;
- hard-code gene candidate thành causal;
- ghi credential/API key;
- thêm remote write không có authorization boundary.

---

## 4. Cấu trúc trả lời sau patch

```yaml
issue_id: GEN-...
status: completed | partial | blocked
baseline_tests:
  passed: ...
final_tests:
  passed: ...
files_changed: []
public_contract_changes: []
invariants_verified: []
known_limitations: []
next_recommended_issue: ...
```

Nếu partial phải nói rõ phần nào chưa đạt acceptance criteria; không dùng “completed” theo cảm giác.

---

## 5. Branch/commit convention đề xuất

```text
branch: gen/<issue-id>-<short-name>
commit: <issue-id>: <imperative summary>
```

Ví dụ:

```text
GEN-021: add two-locus independent assortment plugin
```

---

## 6. Definition of Done chung

- Acceptance criteria của issue đều pass.
- Test mới kiểm failure mode, không chỉ happy path.
- Probability sum 1 trong tolerance.
- Output JSON serializable.
- Error machine-readable.
- Evidence language không nâng.
- Docs/schema đồng bộ.
- Changelog/version đúng.
- No secret/network trong core.

---

## 7. Review checklist khoa học

- Model set có bỏ sót `non_genetic`/UNK không?
- Dữ liệu cha mẹ/đời con có độc lập không?
- Có ascertainment bias không?
- Có genotype-specific survival không?
- Có phụ thuộc tuổi/mùa lông không?
- Reciprocal cross có được giữ riêng không?
- Sex có phương pháp xác định không?
- Câu chữ có vượt G-level không?

---

## 8. Review checklist code

- Circular import?
- Float normalization?
- Zero probability handling?
- Log-space stability?
- Pedigree cycle?
- Unknown founder assumption?
- Same individual train/test leakage nếu ML?
- Hard block evaluated before sorting?
- Input mutation?
- Deterministic output order?
- Version snapshot?

---

## 9. Handoff theo pha

### Model A — Domain/contract

Làm schema, dataclass, validation, migration. Không viết model statistics.

### Model B — Genetics core

Làm genotype/gamete/cross/plugin. Không viết database/UI.

### Model C — Statistical inference

Làm likelihood/Bayes/calibration/simulation study. Không sửa phenotype ontology.

### Model D — Pedigree/population

Làm relationship/inbreeding/diversity. Không sửa cross formula.

### Model E — Image phenotype

Làm calibrated feature pipeline và OOD. Không xuất gene.

### Model F — Platform

Làm API/database/audit/RBAC. Không đưa network vào core.

Agent tổng hợp phải review contract giữa các model; không merge tự động khi schema khác nhau.

---

## 10. First recommended implementation sequence

1. `GEN-001` baseline/golden outputs.
2. `GEN-002` genotype object + property tests.
3. `GEN-010` breeding record schema.
4. `GEN-011` simulation recovery harness.
5. `GEN-020` plugin protocol.
6. `GEN-021` two-locus independent model.
7. `GEN-030` pedigree governance extension.
8. `GEN-040` Pareto mate optimizer.
9. `GEN-050` API/database only after core contracts freeze.

---

## 11. Minimal handoff package to another model

Cần gửi ít nhất:

```text
SKILL.md
references/scientific-boundaries.md
references/architecture.md
references/data-contracts.md
references/algorithms.md
planning/backlog.json
src/
tests/
examples/
```

Không chỉ gửi một đoạn code rời vì sẽ mất evidence boundaries.


<!-- END MODEL CODER HANDOFF -->