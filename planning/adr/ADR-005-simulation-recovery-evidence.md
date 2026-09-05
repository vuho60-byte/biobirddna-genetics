# ADR-005 — Recovery, calibration và misspecification là ba metric riêng

## Status

Accepted.

## Context

Posterior chỉ tương đối trong candidate set. Top-1 accuracy không phải calibration,
và posterior cao không chứng minh candidate set có absolute fit tốt. Parent
configuration còn là tham số rời rạc được chọn từ dữ liệu.

## Decision

- Recovery top-1 yêu cầu unique posterior maximum; tie được báo thành tập, không
  phá tie theo thứ tự input.
- Calibration chỉ dùng ensemble design có nhiều realized true models và so
  posterior với true-model indicator.
- Misspecification dùng multinomial deviance với parametric bootstrap; parent
  configuration được refit trong mỗi bootstrap replicate.
- Candidate set luôn giữ đúng một `M0_NONGENETIC` và ít nhất một model đối chứng
  khác; generator `non_genetic` chỉ nhận genotype `NA`.
- Bootstrap resolution, snapshot consistency, parent phenotype/genotype và work
  budget đều fail closed trước sampling.
- Python và Node dùng cùng SHA-256 seed 53-bit, LCG modulo `2^53` và UTF-8 bucket
  ordering; cùng request/version phải tạo cùng recovery output ở mức parsed-JSON
  với float tolerance (rel_tol `1e-9`, abs_tol `1e-12`). Byte-exact chỉ được bảo
  đảm cho double-run trên cùng một máy — libm khác platform có thể lệch ULP.
  Output giữ thứ tự `candidate_model_ids` đầu vào, còn tie policy độc lập với
  thứ tự đó.
- Mọi output là synthetic decision support, không nâng evidence grade hoặc xác
  nhận gene.

## Consequences

- Recovery, discrimination, calibration và absolute-fit warning có field riêng.
- Cùng seed và contract tạo output deterministic, nhưng không hứa ổn định byte
  qua thay đổi thuật toán/version.
- Cross-runtime equality là invariant trong cùng engine version ở mức parsed-JSON
  float-tolerance (byte-exact chỉ trên cùng máy); đổi PRNG hoặc seed material là
  thay đổi thuật toán cần golden/changelog mới.
- GEN-012 phải tách survival khỏi segregation; GEN-011 không giả lập
  genotype-specific survival.
