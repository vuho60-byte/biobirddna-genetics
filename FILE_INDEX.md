# File index

## Entry points

- `SKILL.md` — operational skill contract.
- `README.md` — quick start and release verification.
- `agents/openai.yaml` — Codex interface metadata.
- `scripts/biobirddna_cli.py` — local CLI.
- `scripts/atomic_rename.py` — platform-specific atomic no-replace finalizer.
- `scripts/validate_package.py` — complete development and release gate.
- `scripts/verify_examples.py` — golden-output and probability verifier.
- `scripts/verify_package.py` — exact inventory, checksum, and version verifier.
- `scripts/rebuild_release_inventory.py` — deterministic checksum inventory builder.
- `npm/package.json` — public scoped release-candidate metadata for the dependency-free ESM package.
- `npm/bin/biobirddna-genetics.mjs` — Node CLI for all six public operations.
- `npm/bin/biobirddna-feedback.mjs` — immutable, proposal-only feedback artifact CLI.

## Architecture and science

- `references/architecture.md`
- `references/algorithms.md`
- `references/data-contracts.md`
- `references/scientific-boundaries.md`
- `references/t02-integration.md`
- `references/implementation-roadmap.md`
- `references/model-card.md`
- `references/autoresearch-feedback.md`
- `references/postgres-schema-draft.sql`
- `references/openapi-draft.yaml`

## Planning for coding models

- `planning/MODEL_CODER_HANDOFF.md`
- `planning/backlog.json`
- `planning/adr/*.md`
- `planning/adr/ADR-007-oauth-mcp-cloud-adapter.md` — stateless OAuth transport boundary.

## Code

- `src/biobirddna_genetics/` — modular reference implementation.
- `npm/src/` — pure Node.js ESM implementation and feedback adapter.
- `npm/test/` — Node parity, safety and packed-install smoke tests; not shipped.
- `standalone/biobirddna_genetics_standalone.py` — independent Python implementation.
- `standalone/biobirddna_genetics_node.mjs` — dependency-free Node.js cross simulator.
- `tests/` — unit, integration, golden, package, installer, and release-inventory tests.

## Contracts and examples

- `assets/models/default-models.json`
- `assets/schemas/`
- `assets/templates/`
- `examples/input/`
- `examples/output/`

## Knowledge dependency

- `knowledge/T02_knowledge_base_v1.0.json`

## Release metadata

- `BioBirdDNA_Genetics_Skill_Blueprint_v0.1.0.md` — original architecture blueprint.
- `TEST_REPORT.md` — verified test scope and results.
- `MANIFEST.json` — runtime/package metadata.
- `PACKAGE_CHECKSUMS.json` — authoritative exact-file SHA-256 manifest.
- `SHA256SUMS.txt` — portable checksum list.
- `CHANGELOG.md`, `VERSION` — release history and version.

The npm package includes `npm/LICENSE` with `UNLICENSED`; no public distribution
license is implied. This tree is a public release candidate, not evidence that a
registry publication has occurred.
The generated `.tgz` is a sibling release artifact and is intentionally not
embedded in the exact skill tree.

## Utility commands

- `run_tests.sh`, `run_tests.bat` — run Python and npm test suites.
- `run_demo.sh`, `run_demo.bat` — generate all example outputs.
- `scripts/install_skill.sh`, `scripts/install_skill.ps1` — verified staged installation.
