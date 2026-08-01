# ADR-001 — Deterministic core không phụ thuộc framework

## Status

Accepted.

## Context

Punnett, genotype parsing và pedigree là logic toán học cần tái lập, audit và chạy được trên nhiều môi trường.

## Decision

Core Python v0.1.0 chỉ dùng standard library và không import network/database/UI.

## Consequences

- Dễ chạy và review.
- Chưa có validation framework mạnh.
- Adapter production được thêm ngoài core.
