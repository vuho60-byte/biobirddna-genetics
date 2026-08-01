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

Prototype v0.2.0 chưa bắt buộc hai field này trong mọi CLI operation để giữ backward compatibility, nhưng production envelope phải thêm.

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
- `normal_label != target_label`.
- với `incomplete_dominant`, `intermediate_label` khác cả `normal_label` và
  `target_label`.
- probability fields trong `[0,1]`.
- `model_id` + model version production phải immutable.
- input phải đúng JSON type; runtime không đổi `null`, boolean hoặc numeric string
  thành số/chuỗi hợp lệ.

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
Các giá trị là trọng số tương đối hữu hạn nên có thể lớn hơn 1. Không được trộn
hai mode known/uncertain trong cùng request. Với model `non_genetic`, known mode
phải dùng `NA`/`NA`, còn uncertain mode chỉ được có key `NA` với trọng số dương.

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

Nếu truyền `models` tường minh, tập này phải có ít nhất hai model và chứa đúng
một model canonical `{"model_id":"M0_NONGENETIC","mode":"non_genetic"}`.
Không cho phép inference một-model vì posterior `1.0` khi đó chỉ là hệ quả của
candidate set, không phải bằng chứng. Mỗi `model_id` phải duy nhất; posterior bằng
nhau được sắp xếp ổn định theo `model_id`.

`model_priors` là trọng số tương đối hữu hạn, không âm; không bắt buộc từng trọng
số nhỏ hơn hoặc bằng 1, nhưng tổng phải lớn hơn 0. JSON input phải là JSON chuẩn:
`NaN`, `Infinity` và `-Infinity` đều bị từ chối.

### Sex mode

Một request phải dùng một trong hai kiểu:

- tất cả `male/female`;
- tất cả `unknown`.

Không trộn known và unknown sex trong một request v0.1.0.

### Counts

- integer trong `[0, 1_000_000]`; đây là work budget chung của Python/Node;
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
Model không có cấu hình bố mẹ tương thích trả `log_evidence:null` và posterior
`0` khi còn ít nhất một model có evidence hữu hạn; runtime không xuất số
`NaN`/`Infinity`.

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
Mọi hypothesis weight là số hữu hạn không âm và tổng phải >0; giá trị có thể lớn
hơn 1 vì đây là trọng số tương đối. Năm pairing weights cũng phải hữu hạn, không
âm. `score` không hữu hạn là lỗi fail-closed, không được serialize thành `null`.
Cùng một `individual_id` ở cả `males[]` và `females[]` luôn thêm
`same_individual_pair` vào `block_reasons`, bất kể
`max_offspring_inbreeding` hoặc `block_direct_ancestor_and_full_sibling`.

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
  "needs_recapture": false,
  "limits": {
    "measurement_domain": "visible_phenotype_only",
    "biochemical_pigment_status": "unconfirmed",
    "gene_status": "not_assessed",
    "inheritance_status": "not_assessed"
  }
}
```

`score` là normalized heuristic, không phải gene posterior.

`image_quality` Q0 hoặc Q1 là hard abstention: runtime phải bỏ qua mọi feature
nhìn thấy dù score đầu vào mạnh, trả `primary_hypothesis:"UNK"`,
`needs_recapture:true`, không đưa giả thuyết thay thế và giữ nguyên toàn bộ giới
hạn phenotype-only. Q2 trở lên trả `needs_recapture:false`; Q-level không nâng
evidence grade và không tạo gene/genotype claim.

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
  "engine_version": "0.2.0",
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
- `assets/schemas/breeding-record.schema.json`
- `assets/schemas/simulation-recovery-request.schema.json`
- `assets/schemas/feedback-input.schema.json`
- `assets/schemas/feedback-artifact.schema.json`

Prototype không phụ thuộc thư viện JSON Schema; schema hiện phục vụ contract, validation ngoài và model handoff.

---

## 19. Normalized breeding-record bundle (GEN-010)

Schema: `assets/schemas/breeding-record.schema.json`.

Contract bắt buộc:

- `pairings[]` giữ rõ `male_parent_id` và `female_parent_id`, nguồn xác định giới
  tính và confidence của từng bố mẹ, cùng `reciprocal_group_id` khi có;
- `clutches[].offspring_records_coverage` là `complete`, `partial` hoặc `unknown`;
- bốn `stage_counts` luôn xuất hiện: `eggs_laid`, `fertile_eggs`, `hatched`,
  `survived_to_first_molt`;
- mỗi stage phân biệt `status:"missing"` với
  `status:"observed", count:0`;
- mỗi hồ sơ trứng/con ghi riêng `fertility_status`, `hatch_status`,
  `survival_status`, `outcome`, sex provenance và phenotype timeline;
- `sequence_index` xác định thứ tự timeline, không dùng tên stage để suy thứ tự;
- khi coverage hồ sơ và stage đều `complete`, runtime phải đối soát count với
  từng hồ sơ và fail closed nếu mâu thuẫn;
- `outcome` phải phù hợp đường fertility → hatch → survival;
- output chỉ là normalized observations, giữ `gene_status:"not_assessed"`.

Không được suy gene/genotype từ `parent_phenotype` hoặc `phenotype_timeline`.

---

## 20. Synthetic simulation-recovery request (GEN-011)

Schema: `assets/schemas/simulation-recovery-request.schema.json`.

Mỗi experiment khai báo snapshot model sinh dữ liệu, genotype và phenotype bố mẹ,
brood sizes, số replicate và candidate-model IDs. Candidate set phải chứa đúng
một `M0_NONGENETIC` và ít nhất một model đối chứng khác. Với generating model
`non_genetic`, cả `male_genotype` và `female_genotype` phải đúng chuỗi `NA`.

Contract output:

- seed là 53 bit từ 64 bit đầu của SHA-256 trên chuỗi
  `base_seed|experiment_id|brood_size|replicate_index`, mask về `2^53-1`;
- Python và Node dùng cùng LCG modulo `2^53`
  (`a=6364136223846793005`, `c=1442695040888963407`) và sort bucket theo byte
  UTF-8, nên cùng request/version tạo cùng stream và output;
- `candidate_model_ids` giữ nguyên thứ tự input trong output; tie không dùng thứ
  tự đó để chọn winner;
- posterior luôn mang nghĩa tương đối trong candidate set;
- `best_models` chứa toàn bộ tie set; `best_model` chỉ có giá trị khi unique;
- `top1_rate` chỉ tính unique top-1, còn `top_set_contains_true_rate` là metric
  tie-inclusive riêng;
- uncertainty gồm Wilson interval, posterior quantiles, entropy và Brier score;
- calibration chỉ được ước lượng trên design có ít nhất hai realized true models
  và so posterior với true-model indicator, không so với top-1 accuracy;
- absolute-fit dùng multinomial deviance và parametric bootstrap, refit parent
  configuration trong từng bootstrap replicate;
- `candidate_model_set_misspecified_suspected` chỉ xuất hiện khi mọi candidate
  có bootstrap p-value dưới alpha;
- `B + 1 > 1 / alpha` và work budgets phải qua trước khi sampling;
- model snapshot, genotype và phenotype bố mẹ phải tương thích nội bộ;
- identifier, description, allele, genotype và phenotype label không được có whitespace ở
  đầu/cuối; runtime không tự trim để biến input sai thành model khác;
- ID model built-in (`M0_NONGENETIC`, `M1_AR`–`M5_ZD`) luôn giữ semantics
  canonical kể cả khi model sinh dữ liệu bị loại khỏi candidate set; một generator
  misspec tùy biến phải dùng ID riêng;
- output luôn có `synthetic_only:true`, `decision_support_only:true`,
  `evidence_grade:"G0"` và `gene_status:"not_assessed"`.

Harness này chưa mô hình genotype-specific survival; phần đó thuộc GEN-012.

---

## 21. OAuth MCP transport contract (GEN-015)

Remote endpoint công khai đúng bảy tên tool:

1. `biobirddna_simulate_cross`
2. `biobirddna_infer_candidate_models`
3. `biobirddna_screen_visible_pigment`
4. `biobirddna_rank_candidate_pairs`
5. `biobirddna_normalize_breeding_records`
6. `biobirddna_simulate_model_recovery`
7. `biobirddna_prepare_feedback`

Mỗi tool là read-only, idempotent, closed-world và non-destructive. Sáu compute
tool nhận đúng schema nguồn tương ứng và dispatch vào fixed operation; không có
tool tên tùy ý, raw code, SQL, filesystem path hoặc generic operation selector.

`biobirddna_prepare_feedback` nhận feedback input không chứa trusted package
integrity. Server gắn integrity đã pin rồi tạo artifact; client không được ghi đè
hash, `review_status`, `auto_apply` hoặc `target_files`. Output bắt buộc giữ
`pending_human_review`, `auto_apply:false`, `target_files:[]` và evidence/limits.

HTTP request chỉ được dispatch sau khi access token qua issuer, signature, expiry,
exact MCP resource audience và allowed-user check. Response/log/evidence không
được chứa password, access/refresh token, biological payload hoặc connection ID
cá nhân. Adapter không tạo table, bucket, queue hay model-run persistence.
