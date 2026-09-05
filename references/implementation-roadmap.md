# Lộ trình triển khai production

Tài liệu này biến kiến trúc thành các pha mà model/coder khác có thể thực hiện tuần tự.

---

## Pha 0 — Khóa invariant và test nền

### Mục tiêu

Bảo đảm mọi model khác hiểu cùng một semantics trước khi thêm tính năng.

### Việc làm

- Đọc `SKILL.md` và `scientific-boundaries.md`.
- Chạy toàn bộ unit test.
- Tạo golden output cho năm example.
- Ghi SHA-256 manifest.
- Không thêm framework/network.

### DoD

- test sạch;
- demo chạy;
- không có claim vượt evidence;
- output sum probability = 1.

---

## Pha 1 — Harden deterministic core

### Scope

- robust genotype object thay string nội bộ;
- exact rational/Fraction option cho Punnett;
- custom sex ratio autosomal;
- penetrance theo genotype/sex;
- phenotype state tách carrier visible status;
- property tests.

### Không làm

- database;
- image ML;
- multi-locus.

### DoD

- backward compatible JSON;
- benchmark 10.000 cross;
- no floating drift > tolerance;
- mutation tests cho Z-linked.

---

## Pha 2 — Breeding records và inference v0.2

### Trạng thái v0.3.0

GEN-010 và GEN-011 đã hoàn thành: normalized breeding records, seeded recovery,
design-level calibration và absolute-fit misspecification guard. GEN-013 và
GEN-014 cũng đã hoàn thành với Node ESM parity, packed-install gate và immutable
feedback artifact chờ người duyệt.
Genotype-specific survival vẫn chưa được triển khai và thuộc GEN-012; feedback
adapter không có quyền gọi hoặc tự apply AutoResearch.

### Scope

- `pairings`, `clutches`, `offspring` schema;
- sex method/confidence;
- multiple broods;
- survival stages;
- missing category semantics;
- configurable model priors;
- Bayes factor/log evidence report;
- parent genotype posterior export.

### DoD

- ingest ít nhất 20 synthetic pedigrees;
- recover known generating model trong simulation study;
- calibration report;
- abstain khi model set misspecified.

---

## Pha 3 — Multi-locus plugin

### Scope

- generic locus object;
- independent assortment;
- two-locus epistasis phenotype map;
- optional linkage/recombination fraction;
- plugin registry;
- JSON model DSL.

### DoD

- reproduce known synthetic ratios;
- reject invalid phenotype maps;
- complexity guard;
- no gene-name claim.

---

## Pha 4 — Pedigree/population governance

### Scope

- unknown-founder uncertainty;
- half-sibling/first-cousin labels;
- founder contribution;
- mean kinship;
- effective population size proxy;
- limit repeated parent contribution;
- genomic relationship adapter future interface.

### DoD

- test pedigrees with loops rejected;
- relationship matrix compared with independent implementation;
- hard block cannot be overridden by weights.

---

## Pha 5 — Mate optimizer v0.2

### Scope

- separate research and production objectives;
- Pareto frontier instead of only scalar score;
- multi-generation simulation;
- diversity constraints;
- capacity constraints per breeder;
- explainable ranking.

### DoD

- output eligible, blocked, Pareto set;
- EIG validated by simulation;
- report sensitivity to weights;
- no auto-execution.

---

## Pha 6 — Image phenotype pipeline

### Scope

- media provenance and consent gate;
- RAW/JPEG QA;
- calibration target detection;
- bird segmentation;
- body region segmentation;
- color feature extraction;
- longitudinal alignment;
- OOD/UNK;
- two-reviewer labeling.

### Dependency candidates

Chọn sau technical review; không khóa ở blueprint. Core phải nhận feature chuẩn, không phụ thuộc CV library cụ thể.

### DoD

- split theo individual/capture session;
- external device validation;
- feature repeatability;
- model card;
- no gene output.

---

## Pha 7 — API/database

### Scope

- Postgres/Supabase schema;
- migration;
- API adapter;
- auth/RBAC;
- idempotency;
- immutable run registry;
- audit event;
- job queue.

### DoD

- contract tests;
- role tests;
- consent withdrawal path;
- deterministic core remains network-free;
- exact run reproducibility.

---

## Pha 8 — Molecular integration

### Scope

- sample/aliquot IDs;
- chain of custody;
- variant import;
- genotype QC;
- cosegregation;
- family-based association;
- genomic kinship;
- G5 evidence report.

### DoD

- sample swap checks;
- missing genotype handling;
- no G5 without lineage/QC;
- publication language guard.

---

## Pha 9 — Governance và v1.0

### Scope

- scientific review board;
- veterinary review;
- privacy/legal review;
- model registry lifecycle;
- incident response;
- external replication;
- public documentation.

### DoD

- validated dataset;
- model cards;
- data sheets;
- audit trail;
- reproducible deployment;
- documented known failure modes.

---

## Ưu tiên ngay sau v0.3.0

1. GEN-012: tách segregation, hatch và genotype-specific survival.
2. GEN-030: mở rộng pedigree governance và founder uncertainty.
3. GEN-040: tách research pairing và target pairing.
4. GEN-041: thêm Pareto ranking và weight sensitivity.
5. Chỉ tạo database/importer sau khi core contracts ổn định.
