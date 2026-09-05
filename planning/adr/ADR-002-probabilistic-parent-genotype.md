# ADR-002 — Kiểu gene cha mẹ được biểu diễn bằng posterior

## Status

Accepted.

## Context

Phần lớn chim không có DNA/test cross đủ để gán genotype chắc chắn.

## Decision

Genotype posterior là first-class input. Engine marginalize qua mọi cấu hình thay vì lấy argmax mặc định.

## Consequences

- Output phản ánh bất định tốt hơn.
- Input phức tạp hơn.
- Optimizer phải tính entropy genotype.
