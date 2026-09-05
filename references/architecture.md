# Kiến trúc sâu — BioBirdDNA Genetics Skill

**Mã kiến trúc:** `CBIM-DNA-GEN-ARCH`  
**Phiên bản:** `0.3.0`  
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
- **Portable:** lõi v0.3.0 chạy bằng Python standard library; npm adapter không có runtime dependency.
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
same_individual_pair
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
- Q0/Q1 image → `UNK` + `needs_recapture:true` dù feature score mạnh.

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

### V0.2 — Breeding records và model recovery

- normalized pairing/clutch/offspring records;
- explicit missing-versus-observed-zero semantics;
- seeded synthetic recovery and design-level calibration;
- absolute-fit misspecification guard with `M0_NONGENETIC`;
- Python and npm/Node adapters under the same evidence boundary, cùng seed
  SHA-256/LCG và exact simulation stream trong một engine version;
- Node feedback adapter chỉ tạo immutable review artifact; không gọi hoặc ghi vào
  AutoResearch, không tái dùng health gate 3.0, fail closed với disguised
  mutation/claim/score và không tự apply.

Xem `planning/adr/ADR-005-simulation-recovery-evidence.md` và
`planning/adr/ADR-006-node-package-feedback-boundary.md`.

Two-locus/epistasis và genotype-specific survival vẫn nằm ở backlog sau v0.3.0.

### V0.2 adapter addendum — OAuth MCP (GEN-015)

- dedicated Supabase project, table-free và stateless;
- đúng bảy MCP tool cố định: sáu compute operation cộng một feedback transformer;
- OAuth 2.1 kiểm issuer, exact resource audience và allowed user trước dispatch;
- Edge chỉ gọi fixed Node operations từ snapshot đã kiểm checksum, không có generic
  dispatch hoặc network trong scientific core;
- feedback luôn `pending_human_review`, `auto_apply:false`, `target_files:[]`;
- Node updater/watcher fail closed, yêu cầu security pin cho sensitive change và
  rollback về immutable snapshot khi live verification thất bại;
- consent UI chạy trên dedicated Cloudflare Worker, không lưu biological payload.

Adapter này không triển khai persistence, stored-record RBAC, consent cho dữ liệu
đã lưu, welfare override hoặc publication flow. Các capability đó vẫn thuộc
GEN-050/GEN-051. Xem `planning/adr/ADR-007-oauth-mcp-cloud-adapter.md`.

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
8. **Không runtime dependency ngoài ở v0.3.0.** Dễ chạy, audit và chuyển giao.

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
