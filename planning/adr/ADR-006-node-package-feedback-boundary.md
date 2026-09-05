# ADR-006 — Node package parity và feedback chỉ dừng ở review artifact

## Status

Accepted.

## Context

BioBirdDNA cần chạy được trong worker Node mà không gọi Python, đồng thời cần một
đường đưa tín hiệu production vào quy trình cải tiến. AutoResearch hiện hữu có
health gate, proposal ledger và quyền promotion riêng; tái sử dụng score 3.0 hoặc
cho adapter tự ghi/apply sẽ trộn evidence sinh học với health vận hành.

## Decision

- Phát hành package ESM private `@biobirddna/genetics`, Node.js >=20, zero
  dependencies và giữ public-operation parity với Python modular.
- Chỉ tạo tarball/cài từ file local trong v0.3.0; không publish registry và giữ
  trạng thái `UNLICENSED`.
- `biobirddna-feedback` nhận normalized production-review envelope: raw operation
  output có thể nằm trong `result.output`, còn `result.evidence` và
  `result.limits` là metadata review bắt buộc.
- Adapter xác minh provenance/integrity, tạo ID SHA-256 deterministic và chỉ xuất
  artifact `pending_human_review`, `auto_apply:false`, `target_files:[]`.
- Adapter không gọi AutoResearch, không ghi proposal ledger, không lấy run-lock,
  không dùng health gate 3.0 làm score BioBirdDNA và không tự approve/apply.
- Mutation intent, molecular confirmation dưới G5, sai version/hash/integrity và
  envelope không có evidence boundary đều fail closed.
- Mutation/claim/health-score guard canonicalize snake, camel và Pascal case;
  structured gene hoặc variant claim có trạng thái confirmed/identified/causal
  dưới G5 và AutoResearch health-score trá hình đều bị chặn.
- Mọi raw public-operation output được đặt trong `result.output` phải là cây JSON
  hữu hạn; trạng thái log-evidence bất khả thi dùng `null`, không dùng
  `NaN`/`Infinity`.

## Consequences

- Worker Node có thể dùng cùng sáu operation mà không cần cây skill ghi được hoặc
  Python runtime.
- Raw output vẫn audit được trong production envelope nhưng không được sao chép
  thành numerical improvement score trong feedback artifact.
- Chuyển artifact thành proposal thật luôn là bước riêng có reviewer người và
  tuân theo gate của AutoResearch hiện hữu.
- Bất kỳ AutoResearch integration có quyền ghi/promotion trong tương lai cần ADR
  và authorization mới; ADR này không cấp quyền đó.
