# Báo cáo hồ sơ cặp ghép

Hồ sơ: SYNTHETIC-BREEDER-001. BioBirdDNA 0.4.0.

Nguồn khai báo: Synthetic example for software validation only.
Dữ liệu minh họa tổng hợp (synthetic), không phải đàn chim thật.

**Chỉ hỗ trợ xem xét; cần người có trách nhiệm duyệt trước khi ghép thật.**

## Bố mẹ và đời con

| Cặp | Bố: mã / kiểu hình | Mẹ: mã / kiểu hình | Nguồn giới tính bố / mẹ |
|---|---|---|---|
| CM-MAT-DEMO-010-A | CM-IND-A / normal | CM-IND-B / target\_colour | dna / dna |
| CM-MAT-DEMO-010-B | CM-IND-A / normal | CM-IND-C / normal | dna / dna |

### Lứa CM-CLU-DEMO-010 — cặp CM-MAT-DEMO-010-A

Độ đầy đủ hồ sơ con: complete.

- Trứng đẻ: 3; độ đầy đủ: complete.
- Trứng có phôi: 2; độ đầy đủ: complete.
- Đã nở: 1; độ đầy đủ: complete.
- Sống tới thay lông đầu: chưa ghi nhận; độ đầy đủ: unknown.

### Lứa SYNTHETIC-EMPTY-CLUTCH — cặp CM-MAT-DEMO-010-B

Độ đầy đủ hồ sơ con: complete.

- Trứng đẻ: 0; độ đầy đủ: complete.
- Trứng có phôi: 0; độ đầy đủ: complete.
- Đã nở: 0; độ đầy đủ: complete.
- Sống tới thay lông đầu: 0; độ đầy đủ: complete.

| Mã trứng/con | Lứa | Kết quả | Giới tính | Diễn biến kiểu hình |
|---|---|---|---|---|
| CM-OFF-DEMO-010-A | CM-CLU-DEMO-010 | Mất theo dõi | chưa rõ | 0 / first\_plumage: normal |
| CM-OFF-DEMO-010-B | CM-CLU-DEMO-010 | Không thụ tinh | chưa rõ | chưa ghi nhận |
| CM-OFF-DEMO-010-C | CM-CLU-DEMO-010 | Chết phôi | chưa rõ | chưa ghi nhận |

## Mô phỏng và xác suất kiểu hình

Tên màu và số con quan sát không tự cho biết kiểu gene. Chỉ tính khi có mô hình và kiểu gene giả định được nhập riêng.
Xác suất dưới đây tính trên toàn bộ con, gộp cả giới tính (không phải tỷ lệ có điều kiện trong riêng trống hoặc mái).
Ví dụ 25% mái màu mục tiêu nghĩa là 25% tổng số con; một lứa nhỏ không được bảo đảm đúng tỷ lệ này.

### Cặp CM-MAT-DEMO-010-A — mô hình giả định M1\_AR

Giả định đầu vào (allele ký hiệu, không phải tên gene đã xác nhận):

```json
{
  "model": {
    "model_id": "M1_AR",
    "mode": "autosomal_recessive",
    "normal_allele": "N",
    "variant_allele": "m",
    "normal_label": "normal",
    "target_label": "target_colour",
    "penetrance": 1.0,
    "observation_error": 0.01
  },
  "male_genotype": "N/m",
  "female_genotype": "m/m"
}
```

| Giới tính và kiểu hình | Xác suất trên tổng số con |
|---|---:|
| mái: normal | 25.00% |
| mái: target\_colour | 25.00% |
| trống: normal | 25.00% |
| trống: target\_colour | 25.00% |

## Cận huyết và khuyến nghị xem xét cặp

F tính từ phả hệ khai báo. Thiếu tổ tiên không đồng nghĩa không có quan hệ; F = 0 không chứng minh không cận huyết.
Ngưỡng chặn mặc định của engine là F > 6,25%; đây không phải ngưỡng an toàn sức khỏe được chứng nhận.
Nếu chưa có cờ phúc lợi, sức khỏe vẫn là chưa đánh giá. Ưu tiên hoàn thiện hồ sơ và để người phụ trách/thú y xem xét.

| Cặp | F kỳ vọng của con | Khuyến nghị |
|---|---:|---|
| CM-MAT-DEMO-010-A | 25.00% | DỪNG xem xét ghép: Anh chị em ruột; Hệ số F vượt ngưỡng chặn mặc định |
| CM-MAT-DEMO-010-B | 0.00% | Chưa phê duyệt ghép; bổ sung/phúc tra phả hệ, sức khỏe và chờ người phụ trách duyệt |

## Dữ liệu thiếu và giới hạn bằng chứng

Evidence grade: `not_assessed`. Không tự cấp G-level. Posterior/mô phỏng không xác nhận gene.
Số 0 đã quan sát khác dữ liệu chưa ghi nhận. Giữ cả trứng không nở, chết phôi, chết sau nở và mất theo dõi.
Mô phỏng chưa mô hình hóa sống sót theo genotype, cấu trúc quần thể hoặc nhiều locus.

- clutches/CM-CLU-DEMO-010/stage\_counts/survived\_to\_first\_molt: survival\_stage\_not\_observed.

### Toàn bộ Limits từ kết quả core và báo cáo

```json
{
  "records": {
    "purpose": "breeding_record_normalization_only",
    "gene_status": "not_assessed",
    "inheritance_status": "not_assessed",
    "decision_support_only": true,
    "sex_method_recorded_not_verified": true,
    "claim_language": "observation_record"
  },
  "simulations": [
    {
      "gene_status": "model_assumption_only",
      "inheritance_status": "simulation_not_confirmation",
      "claim_language": "simulate_cross"
    }
  ],
  "report": {
    "decision_support_only": true,
    "gene_status": "not_assessed",
    "human_review_status": "pending",
    "pedigree_declared_not_verified": true,
    "unknown_founders_assumed_unrelated_for_calculation_only": true,
    "max_offspring_inbreeding": 0.0625,
    "health_safety_not_assessed": true
  }
}
```

### Hồ sơ chuẩn hóa để đối chiếu (giữ nguyên kết quả core)

```json
{
  "schema_version": "1.0.0",
  "case_id": "SYNTHETIC-BREEDER-001",
  "pairings": [
    {
      "pairing_id": "CM-MAT-DEMO-010-A",
      "male_parent_id": "CM-IND-A",
      "female_parent_id": "CM-IND-B",
      "male_parent_phenotype": "normal",
      "female_parent_phenotype": "target_colour",
      "male_parent_sex_method": "dna",
      "male_parent_sex_confidence": 0.99,
      "female_parent_sex_method": "dna",
      "female_parent_sex_confidence": 0.99,
      "reciprocal_group_id": null,
      "cross_direction_key": "CM-IND-A>CM-IND-B"
    },
    {
      "pairing_id": "CM-MAT-DEMO-010-B",
      "male_parent_id": "CM-IND-A",
      "female_parent_id": "CM-IND-C",
      "male_parent_phenotype": "normal",
      "female_parent_phenotype": "normal",
      "male_parent_sex_method": "dna",
      "male_parent_sex_confidence": 0.99,
      "female_parent_sex_method": "dna",
      "female_parent_sex_confidence": 0.99,
      "reciprocal_group_id": null,
      "cross_direction_key": "CM-IND-A>CM-IND-C"
    }
  ],
  "clutches": [
    {
      "clutch_id": "CM-CLU-DEMO-010",
      "pairing_id": "CM-MAT-DEMO-010-A",
      "offspring_records_coverage": "complete",
      "stage_counts": {
        "eggs_laid": {
          "status": "observed",
          "coverage": "complete",
          "count": 3
        },
        "fertile_eggs": {
          "status": "observed",
          "coverage": "complete",
          "count": 2
        },
        "hatched": {
          "status": "observed",
          "coverage": "complete",
          "count": 1
        },
        "survived_to_first_molt": {
          "status": "missing",
          "coverage": "unknown",
          "count": null
        }
      }
    },
    {
      "clutch_id": "SYNTHETIC-EMPTY-CLUTCH",
      "pairing_id": "CM-MAT-DEMO-010-B",
      "offspring_records_coverage": "complete",
      "stage_counts": {
        "eggs_laid": {
          "status": "observed",
          "coverage": "complete",
          "count": 0
        },
        "fertile_eggs": {
          "status": "observed",
          "coverage": "complete",
          "count": 0
        },
        "hatched": {
          "status": "observed",
          "coverage": "complete",
          "count": 0
        },
        "survived_to_first_molt": {
          "status": "observed",
          "coverage": "complete",
          "count": 0
        }
      }
    }
  ],
  "offspring": [
    {
      "offspring_record_id": "CM-OFF-DEMO-010-A",
      "individual_id": "CM-IND-DEMO-010-A",
      "clutch_id": "CM-CLU-DEMO-010",
      "outcome": "lost_to_followup",
      "fertility_status": "fertile",
      "hatch_status": "hatched",
      "survival_status": "unknown",
      "sex": "unknown",
      "sex_method": "unknown",
      "sex_confidence": null,
      "phenotype_timeline": [
        {
          "observation_id": "CM-PHENO-DEMO-010-A",
          "sequence_index": 0,
          "stage": "first_plumage",
          "status": "observed",
          "phenotype": "normal"
        }
      ]
    },
    {
      "offspring_record_id": "CM-OFF-DEMO-010-B",
      "individual_id": null,
      "clutch_id": "CM-CLU-DEMO-010",
      "outcome": "unfertilized",
      "fertility_status": "unfertilized",
      "hatch_status": "not_applicable",
      "survival_status": "not_applicable",
      "sex": "unknown",
      "sex_method": "unknown",
      "sex_confidence": null,
      "phenotype_timeline": []
    },
    {
      "offspring_record_id": "CM-OFF-DEMO-010-C",
      "individual_id": null,
      "clutch_id": "CM-CLU-DEMO-010",
      "outcome": "embryo_death",
      "fertility_status": "fertile",
      "hatch_status": "not_hatched",
      "survival_status": "not_applicable",
      "sex": "unknown",
      "sex_method": "unknown",
      "sex_confidence": null,
      "phenotype_timeline": []
    }
  ],
  "missing_data": [
    {
      "field": "clutches/CM-CLU-DEMO-010/stage_counts/survived_to_first_molt",
      "impact": "survival_stage_not_observed"
    }
  ],
  "evidence": {
    "source": "submitted_breeding_records",
    "genetics_grade": "not_assessed",
    "record_scope": "normalized_observations_only"
  },
  "limits": {
    "purpose": "breeding_record_normalization_only",
    "gene_status": "not_assessed",
    "inheritance_status": "not_assessed",
    "decision_support_only": true,
    "sex_method_recorded_not_verified": true,
    "claim_language": "observation_record"
  }
}
```

Mã SHA-256 input JSON chuẩn hóa: `f7203159bfc3cf58deee92801a46984534074de5df4f3016aa4b8fd53f03a5f4`.
