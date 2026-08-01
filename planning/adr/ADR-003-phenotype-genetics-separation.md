# ADR-003 — Tách phenotype engine khỏi genetics engine

## Status

Accepted.

## Context

Ảnh chỉ hỗ trợ mô tả màu/pattern; dùng ảnh để gán gene tạo circular evidence và overclaim.

## Decision

Module pigment chỉ xuất T02 working labels và limits. Không tạo genotype posterior hoặc inheritance mode.

## Consequences

- Cần handoff qua breeding records.
- Giảm rủi ro image-to-gene shortcut.
