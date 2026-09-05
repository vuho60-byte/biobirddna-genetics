# Model card — BioBirdDNA Genetics v0.3.0

## Intended use

- Educational/research cross simulation.
- Candidate inheritance-model ranking from complete brood data.
- Pedigree and inbreeding decision support.
- Research-pair ranking.
- T02 phenotype screening from precomputed features.
- Normalization and reconciliation of breeding records.
- Seeded synthetic recovery, design calibration and model-set misspecification checks.

## Out-of-scope use

- Confirming a gene from an image.
- Veterinary diagnosis.
- Automatic breeding order.
- Commercial guarantee of colour outcomes.
- Invasive sample collection guidance.
- Population-wide frequency estimation from social media images.

## Implemented models

- One-locus autosomal recessive.
- One-locus autosomal dominant.
- One-locus incomplete dominance.
- One-locus avian Z-linked recessive.
- One-locus avian Z-linked dominant.
- Simple non-genetic baseline.

## Key assumptions

- Two symbolic alleles.
- Equal transmission.
- No linkage for one-locus models.
- No genotype-specific survival in v0.3.0.
- Simple penetrance.
- Simple observation-error smoothing.
- Unknown founders treated as unrelated in pedigree matrix; this can underestimate true relatedness.

## Known limitations

- Posterior depends on candidate model set and priors.
- Small brood sizes produce unstable ranking.
- Phenotype labels can be age/molt dependent.
- Sex misclassification strongly affects Z-linked inference.
- Non-genetic baseline is intentionally simple.
- Pair score is a configurable decision-support function, not a welfare standard.

## Required human review

- Scientific review for inheritance claims.
- Veterinary/welfare review for real pairings.
- Legal/consent review for data and samples.
- Lab review for molecular claims.

## Evidence language

- G0–G2: exploratory/candidate.
- G3–G4: supported/strong model, no confirmed locus.
- G5: locus/variant supported by cosegregation.
- G6: independently replicated or functionally supported mechanism.
