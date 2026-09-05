import assert from "node:assert/strict";
import test from "node:test";

import { dispatch } from "../src/index.mjs";

// Inputs are built via JSON.parse (the real CLI path) so that keys like
// "__proto__" arrive as own data properties, exactly as JSON.parse creates
// them, instead of being swallowed by object-literal semantics.
const DUNDER_LABELS = ["__proto__", "constructor", "prototype"];
const CONTROL_LABEL = "zzcontrolzz";
const PROTOTYPE_SNAPSHOT = Object.getOwnPropertyNames(Object.prototype).sort();

const inferRequest = (label) => JSON.parse(`{
  "target_label": ${JSON.stringify(label)},
  "male_parent_phenotype": "normal",
  "female_parent_phenotype": "normal",
  "offspring": [
    {"sex": "unknown", "phenotype": "normal", "count": 6},
    {"sex": "unknown", "phenotype": ${JSON.stringify(label)}, "count": 2}
  ]
}`);

const renameLabel = (result, from, to) => JSON.parse(JSON.stringify(result).replaceAll(from, to));

test("infer treats __proto__/constructor/prototype target labels like ordinary labels", () => {
  const control = dispatch("infer", inferRequest(CONTROL_LABEL));
  assert.equal(control.best_model, "M4_ZR");
  const controlJson = JSON.parse(JSON.stringify(control));
  for (const label of DUNDER_LABELS) {
    const result = dispatch("infer", inferRequest(label));
    assert.equal(result.best_model, "M4_ZR", `label ${label} silently corrupted inference`);
    assert.equal(result.best_model_posterior, control.best_model_posterior);
    assert.deepEqual(renameLabel(result, label, CONTROL_LABEL), controlJson);
  }
});

test("cross handles dunder allele names without losing probability mass", () => {
  const request = JSON.parse(`{
    "model": {"model_id": "X", "mode": "autosomal_recessive", "variant_allele": "__proto__", "target_label": "target"},
    "male_genotype": "N/__proto__",
    "female_genotype": "__proto__/__proto__"
  }`);
  const result = dispatch("cross", request);
  assert.deepEqual({ ...result.genotype_probabilities }, {
    "female:N/__proto__": 0.25,
    "female:__proto__/__proto__": 0.25,
    "male:N/__proto__": 0.25,
    "male:__proto__/__proto__": 0.25
  });
  assert.deepEqual({ ...result.phenotype_probabilities }, {
    "female:normal": 0.25,
    "female:target": 0.25,
    "male:normal": 0.25,
    "male:target": 0.25
  });
});

test("model prior lookup never reads inherited Object.prototype properties", () => {
  // "constructor" is absent from model_priors: Python treats it as prior 0;
  // Node must match instead of reading Object.prototype.constructor.
  const request = JSON.parse(`{
    "target_label": "target",
    "male_parent_phenotype": "normal",
    "female_parent_phenotype": "normal",
    "models": [
      {"model_id": "constructor", "mode": "autosomal_recessive", "target_label": "target"},
      {"model_id": "M1_AR", "mode": "autosomal_recessive", "target_label": "target"},
      {"model_id": "M0_NONGENETIC", "mode": "non_genetic", "target_label": "target"}
    ],
    "model_priors": {"M1_AR": 1, "M0_NONGENETIC": 1},
    "offspring": [
      {"sex": "unknown", "phenotype": "normal", "count": 6},
      {"sex": "unknown", "phenotype": "target", "count": 2}
    ]
  }`);
  const result = dispatch("infer", request);
  assert.equal(result.best_model, "M1_AR");
  const zeroed = result.ranking.find((item) => item.model_id === "constructor");
  assert.equal(zeroed.prior, 0);
  assert.equal(zeroed.posterior, 0);
});

test("rank-mates accepts __proto__ as a hypothesis model id", () => {
  const request = JSON.parse(`{
    "pedigree": [ {"individual_id": "A", "sex": "male"}, {"individual_id": "B", "sex": "female"} ],
    "males": [ {"individual_id": "A", "sex": "male", "genotype_posteriors": {"__proto__": {"N/m": 1}}} ],
    "females": [ {"individual_id": "B", "sex": "female", "genotype_posteriors": {"__proto__": {"m/m": 1}}} ],
    "hypotheses": [ {"model": {"model_id": "__proto__", "mode": "autosomal_recessive", "target_label": "target"}, "weight": 1} ]
  }`);
  const pair = dispatch("rank-mates", request).pair_ranking[0];
  assert.equal(pair.eligible, true);
  assert.deepEqual(pair.block_reasons, []);
  assert.equal(pair.score, 0.5);
  assert.equal(pair.target_probability, 0.5);
});

test("Object.prototype is not polluted after processing dunder inputs", () => {
  for (const label of DUNDER_LABELS) dispatch("infer", inferRequest(label));
  dispatch("cross", JSON.parse(`{
    "model": {"model_id": "X", "mode": "autosomal_recessive", "variant_allele": "__proto__", "target_label": "target"},
    "male_genotypes": {"N/__proto__": 1},
    "female_genotypes": {"__proto__/__proto__": 1}
  }`));
  assert.deepEqual(Object.getOwnPropertyNames(Object.prototype).sort(), PROTOTYPE_SNAPSHOT);
  assert.equal(Object.getPrototypeOf({}), Object.prototype);
  const fresh = {};
  assert.equal(fresh.normal, undefined);
  assert.equal(fresh["male:normal"], undefined);
  assert.equal(fresh[CONTROL_LABEL], undefined);
});
