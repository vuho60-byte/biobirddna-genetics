import assert from "node:assert/strict";
import test from "node:test";

import { dispatch } from "../src/index.mjs";

const ar = {
  model_id: "M1_AR",
  mode: "autosomal_recessive",
  target_label: "target"
};

test("model phenotype labels cannot collide", () => {
  assert.throws(
    () => dispatch("cross", {
      model: { ...ar, normal_label: "same", target_label: "same" },
      male_genotype: "N/m",
      female_genotype: "m/m"
    }),
    /phenotype_labels_must_be_distinct/
  );
  assert.throws(
    () => dispatch("cross", {
      model: {
        model_id: "M3_INC",
        mode: "incomplete_dominant",
        normal_label: "normal",
        target_label: "target",
        intermediate_label: "target"
      },
      male_genotype: "N/m",
      female_genotype: "m/m"
    }),
    /phenotype_labels_must_be_distinct/
  );
});

test("uncertain cross rejects null and coercible weights", () => {
  for (const invalid of [null, true, "1"]) {
    assert.throws(
      () => dispatch("cross", {
        model: ar,
        male_genotypes: { "N/N": invalid, "N/m": 1 },
        female_genotypes: { "m/m": 1 }
      }),
      /probability_distribution_invalid/
    );
  }
});

test("dispatch rejects non-finite numbers and unknown cross fields", () => {
  const base = {
    model: ar,
    male_genotype: "N/m",
    female_genotype: "m/m"
  };
  assert.throws(
    () => dispatch("cross", { ...base, unused: Infinity }),
    /input_json_non_finite:number/
  );
  assert.throws(
    () => dispatch("cross", { ...base, unknown: true }),
    /cross_request_unknown_field:unknown/
  );
});

test("rank-mates accepts nonnegative relative hypothesis weights", () => {
  const result = dispatch("rank-mates", {
    pedigree: [
      { individual_id: "M", sex: "male" },
      { individual_id: "F", sex: "female" }
    ],
    males: [{
      individual_id: "M",
      sex: "male",
      genotype_posteriors: { M1_AR: { "N/m": 1 } }
    }],
    females: [{
      individual_id: "F",
      sex: "female",
      genotype_posteriors: { M1_AR: { "m/m": 1 } }
    }],
    hypotheses: [{ model: ar, weight: 10 }]
  });
  assert.equal(result.pair_ranking.length, 1);
  assert.equal(result.pair_ranking[0].eligible, true);
});

test("rank-mates rejects unsafe weights and non-finite scores", () => {
  const base = {
    pedigree: [
      { individual_id: "M", sex: "male" },
      { individual_id: "F", sex: "female" }
    ],
    males: [{
      individual_id: "M",
      sex: "male",
      genotype_posteriors: {
        AR: { "m/m": 1 },
        AD: { "N/N": 1 }
      }
    }],
    females: [{
      individual_id: "F",
      sex: "female",
      genotype_posteriors: {
        AR: { "m/m": 1 },
        AD: { "N/N": 1 }
      }
    }],
    hypotheses: [
      {
        model: {
          model_id: "AR",
          mode: "autosomal_recessive",
          target_label: "target"
        },
        weight: 10
      },
      {
        model: {
          model_id: "AD",
          mode: "autosomal_dominant",
          target_label: "target"
        },
        weight: 10
      }
    ]
  };
  assert.throws(
    () => dispatch("rank-mates", {
      ...base,
      weights: { health: -1 }
    }),
    /health_weight_invalid/
  );
  assert.throws(
    () => dispatch("rank-mates", {
      ...base,
      weights: {
        target: 1.7e308,
        information: 1.7e308,
        inbreeding: 0,
        health: 0,
        uncertainty: 0
      }
    }),
    /pair_score_non_finite:M:F/
  );
});

test("inference count is capped by its operational work budget", () => {
  assert.throws(
    () => dispatch("infer", {
      male_parent_phenotype: "normal",
      female_parent_phenotype: "normal",
      offspring: [{
        sex: "unknown",
        phenotype: "normal",
        count: 1_000_001
      }]
    }),
    /offspring_count_invalid/
  );
});

test("cross model rejects schema-invalid coercions", () => {
  const base = {
    model_id: "M1_AR",
    mode: "autosomal_recessive",
    target_label: "target"
  };
  for (const [change, error] of [
    [{ target_label: null }, /target_label_required/],
    [{ model_id: null }, /model_id_required/],
    [{ penetrance: "0.5" }, /penetrance_invalid/]
  ]) {
    assert.throws(
      () => dispatch("cross", {
        model: { ...base, ...change },
        male_genotype: "N/m",
        female_genotype: "m/m"
      }),
      error
    );
  }
  const { model_id: _omitted, ...missingId } = base;
  assert.throws(
    () => dispatch("cross", {
      model: missingId,
      male_genotype: "N/m",
      female_genotype: "m/m"
    }),
    /model_id_required/
  );
  assert.throws(
    () => dispatch("cross", {
      model: { ...base, unknown: true },
      male_genotype: "N/m",
      female_genotype: "m/m"
    }),
    /model_unknown_field:unknown/
  );
});

test("rank-mates rejects null and coercible controls", () => {
  const base = {
    pedigree: [
      { individual_id: "M", sex: "male" },
      { individual_id: "F", sex: "female" }
    ],
    males: [{
      individual_id: "M",
      sex: "male",
      genotype_posteriors: { M1_AR: { "N/m": 1 } }
    }],
    females: [{
      individual_id: "F",
      sex: "female",
      genotype_posteriors: { M1_AR: { "m/m": 1 } }
    }],
    hypotheses: [{ model: ar, weight: 1 }]
  };
  for (const [change, error] of [
    [{ max_offspring_inbreeding: null }, /max_offspring_inbreeding_invalid/],
    [{ block_direct_ancestor_and_full_sibling: null }, /block_direct_ancestor_and_full_sibling_invalid/],
    [{ hypotheses: [{ model: ar, weight: "1" }] }, /model_weight_invalid/]
  ]) assert.throws(() => dispatch("rank-mates", { ...base, ...change }), error);
});

test("rank-mates rejects invalid profiles and unknown fields", () => {
  const base = {
    pedigree: [
      { individual_id: "M", sex: "male" },
      { individual_id: "F", sex: "female" }
    ],
    males: [{
      individual_id: "M",
      sex: "male",
      genotype_posteriors: { M1_AR: { "N/m": 1 } },
      health_risk: 0,
      welfare_blocked: false
    }],
    females: [{
      individual_id: "F",
      sex: "female",
      genotype_posteriors: { M1_AR: { "m/m": 1 } }
    }],
    hypotheses: [{ model: ar, weight: 1 }],
    weights: {}
  };
  const cases = [
    [{ ...base, males: [{ ...base.males[0], individual_id: null }] }, /individual_id_required/],
    [{ ...base, males: [{ ...base.males[0], health_risk: "0.1" }] }, /health_risk_invalid/],
    [{ ...base, males: [{ ...base.males[0], welfare_blocked: null }] }, /welfare_blocked_invalid/],
    [{ ...base, unknown: true }, /rank_mates_request_unknown_field:unknown/],
    [{ ...base, males: [{ ...base.males[0], unknown: true }] }, /genetic_profile_unknown_field:unknown/],
    [{ ...base, weights: { unknown: 1 } }, /weights_unknown_field:unknown/],
    [{ ...base, hypotheses: [{ ...base.hypotheses[0], unknown: true }] }, /model_hypothesis_unknown_field:unknown/]
  ];
  for (const [request, error] of cases) assert.throws(() => dispatch("rank-mates", request), error);
});
