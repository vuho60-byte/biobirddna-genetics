# @biobirddna/genetics

Dependency-free Node 20+ research primitives for evidence-bounded avian inheritance simulations. Outputs are decision support only: phenotype observations do not infer genes, posterior values are relative to candidate models, and no operation creates a breeding order.

This v0.3.0 package is private and `UNLICENSED`. It is distributed locally;
`private: true` prevents registry publication. No public publishConfig is set.

```sh
npm install --global ./biobirddna-genetics-0.3.0.tgz
```

Use the `biobirddna-genetics` bin with JSON input, or import `dispatch` from the package. Supported operations are `cross`, `infer`, `pigment`, `rank-mates`, `normalize-records`, and `simulate-recovery`.

Inference requires at least two candidate models and exactly one canonical
`M0_NONGENETIC`. Simulation recovery applies the same comparison boundary and
requires `NA` parent genotypes for a non-genetic generator. Public results are
strict JSON values; impossible log-evidence is `null`, never `Infinity`.
Inference observations are capped at `1_000_000` per row. Model phenotype labels
must be distinct; hypothesis and pairing weights must be finite and nonnegative.
Python and Node share the v0.3.0 SHA-256 seed, 53-bit LCG, and UTF-8 sampling
order, so the same simulation request produces the same output.

Pigment screening is phenotype-only. `image_quality` Q0 or Q1 always hard
abstains with `primary_hypothesis: "UNK"` and `needs_recapture: true`, even when
visible feature scores are strong; it never infers a gene or genotype.

`biobirddna-feedback --input <completed-envelope.json> --expected-integrity <sha256> [--output <artifact.json>]` creates an immutable, review-only `biobirddna_feedback_v1` artifact. It is local-only, never creates or applies a proposal, and fails closed if production provenance or evidence boundaries are invalid. Import `createFeedbackArtifact(envelope, expectedIntegrity)` for the same deterministic adapter.

The guard canonicalizes snake/camel/Pascal-case keys and rejects disguised
mutation intent, confirmed/causal/identified gene or variant claims below G5,
and AutoResearch health-score fields. Ordinary domain scores and fractions
remain valid provenance.

The feedback input is a normalized production-review wrapper, not a raw operation result: supply `result.evidence`, `result.limits.decision_support_only: true`, and `result.limits.gene_status: "not_assessed"`. If provenance needs the original public operation payload, put it under `result.output`; the adapter does not copy that payload into the feedback artifact.
