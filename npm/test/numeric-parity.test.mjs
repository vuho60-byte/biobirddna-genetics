// Order independence and determinism of the shared numeric primitives.
// The cross-runtime bit-identity sweep lives in tests/test_npm_parity.py
// (NumericPrimitiveParityTests); this file guards the Node side alone:
// reordering summands, candidate models, or simulation candidates must not
// change a single output bit.
import test from "node:test";
import assert from "node:assert/strict";
import { dispatch, __numericPrimitivesForTest as P } from "../src/index.mjs";

test("logSumExp, entropy and orderIndependentMean are permutation invariant", () => {
  const values = [-1.2345, -0.5, -3.75, -2.25, -0.5000000001, -7.125, -0.4999999999];
  const permutations = [
    [...values].reverse(),
    [values[3], values[0], values[6], values[1], values[5], values[2], values[4]],
  ];
  for (const permuted of permutations) {
    assert.strictEqual(P.logSumExp(permuted), P.logSumExp(values));
    assert.strictEqual(P.orderIndependentMean(permuted), P.orderIndependentMean(values));
    assert.strictEqual(P.neumaierSum([...permuted].sort((a, b) => a - b)), P.neumaierSum([...values].sort((a, b) => a - b)));
  }
  const posterior = { a: 0.527, b: 0.265, c: 0.1, d: 0.05, e: 0.033, f: 0.025 };
  const reordered = Object.fromEntries(Object.entries(posterior).reverse());
  assert.strictEqual(P.entropy(reordered), P.entropy(posterior));
});

test("infer result is bit-identical when candidate model order is reversed", () => {
  const shared = { target_label: "target" };
  const models = [
    { model_id: "M1_AR", mode: "autosomal_recessive", ...shared },
    { model_id: "M2_AD", mode: "autosomal_dominant", ...shared },
    { model_id: "M3_INC", mode: "incomplete_dominant", ...shared },
    { model_id: "M4_ZR", mode: "z_linked_recessive", ...shared },
    { model_id: "M5_ZD", mode: "z_linked_dominant", ...shared },
    { model_id: "M0_NONGENETIC", mode: "non_genetic", ...shared },
  ];
  const base = {
    target_label: "target",
    male_parent_phenotype: "normal",
    female_parent_phenotype: "normal",
    offspring: [
      { sex: "unknown", phenotype: "normal", count: 6 },
      { sex: "unknown", phenotype: "target", count: 2 },
    ],
  };
  const forward = dispatch("infer", { ...base, models });
  const reversed = dispatch("infer", { ...base, models: [...models].reverse() });
  assert.deepStrictEqual(reversed, forward);
});

test("simulate-recovery summaries are identical when tied candidates are reordered", () => {
  const experiment = {
    experiment_id: "SIM-TIE-NODE",
    generating_model: { model_id: "M2_AD", mode: "autosomal_dominant", target_label: "target" },
    male_genotype: "N/N",
    female_genotype: "N/N",
    male_parent_phenotype: "normal",
    female_parent_phenotype: "normal",
    brood_sizes: [16],
    replicates: 2,
    candidate_model_ids: ["M1_AR", "M2_AD", "M3_INC", "M4_ZR", "M5_ZD", "M0_NONGENETIC"],
  };
  const reordered = { ...experiment, candidate_model_ids: ["M3_INC", "M5_ZD", "M2_AD", "M1_AR", "M4_ZR", "M0_NONGENETIC"] };
  const request = (item) => ({ schema_version: "1.0.0", base_seed: 11, goodness_of_fit_replicates: 10, misspecification_alpha: 0.2, experiments: [item] });
  const first = dispatch("simulate-recovery", request(experiment)).experiments[0].by_brood_size[0];
  const second = dispatch("simulate-recovery", request(reordered)).experiments[0].by_brood_size[0];
  assert.deepStrictEqual(second.recovery, first.recovery);
  assert.deepStrictEqual(second.uncertainty, first.uncertainty);
  assert.deepStrictEqual(
    second.replicate_results.map((replicate) => replicate.best_models),
    first.replicate_results.map((replicate) => replicate.best_models),
  );
  assert.deepStrictEqual(
    second.replicate_results.map((replicate) => replicate.model_posteriors),
    first.replicate_results.map((replicate) => replicate.model_posteriors),
  );
});
