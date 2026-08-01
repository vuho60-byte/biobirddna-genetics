# Tích hợp T02 vào BioBirdDNA Genetics Skill

Nguồn máy đọc được bundle tại:

```text
knowledge/T02_knowledge_base_v1.0.json
```

T02 định nghĩa 14 lớp ontology, thang E/C/Q, trường phenotype, timeline và nguyên tắc không suy genotype từ ảnh.

---

## 1. Mapping lớp T02

Các nhãn phenotype screen v0.2.0:

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

- Q0–Q1: không measurement-ready; bắt buộc `UNK` và `needs_recapture:true`,
  không cho feature score vượt cổng chất lượng.
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
