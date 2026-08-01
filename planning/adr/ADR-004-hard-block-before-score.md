# ADR-004 — Hard constraints trước hàm điểm ghép

## Status

Accepted.

## Context

Một scalar score có thể vô tình cho cặp quan hệ gần thắng nhờ xác suất màu cao.

## Decision

Welfare, direct ancestor, full sibling và ngưỡng cận huyết cấu hình được đánh giá trước score. Blocked pair không eligible dù score số cao.

## Consequences

- Report phải giữ blocked list và lý do.
- Weight không thể vô hiệu hóa policy.
