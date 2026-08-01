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

## 10. Recommended implementation sequence

Hoàn thành trong v0.2.0:

1. `GEN-001` baseline/golden outputs.
2. `GEN-010` normalized breeding records.
3. `GEN-011` simulation recovery/calibration harness.
4. `GEN-013` — package Node ESM parity cho sáu operation và local packed-install gate.
5. `GEN-014` — feedback adapter proposal-only theo ADR-006; không có quyền AutoResearch
   mutation/promotion.
6. `GEN-015` — OAuth-protected stateless MCP adapter theo ADR-007; đúng bảy tool,
   exact Edge parity, no persistence và feedback vẫn pending human review.

Thứ tự khuyến nghị tiếp theo phải được đối chiếu lại với dependency trong
`planning/backlog.json`; không tự bypass dependency. Sau GEN-015, ready là
`GEN-002`, `GEN-012`, `GEN-030` và `GEN-060`; các mục còn lại chỉ mở khi dependency
của chính chúng hoàn thành:

1. `GEN-002` genotype object + property tests.
2. `GEN-012` genotype-specific survival/observation stages.
3. `GEN-030` pedigree governance extension.
4. `GEN-020` plugin protocol sau `GEN-002`.
5. `GEN-040` research/production objectives sau `GEN-030`.
6. `GEN-050` API/database chỉ sau khi core contracts ổn định.

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
