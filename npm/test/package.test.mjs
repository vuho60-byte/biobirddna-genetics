import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";
import test from "node:test";

import { createFeedbackArtifact, dispatch, BioBirdDNAError, __sampleForTest } from "../src/index.mjs";

const execFileAsync = promisify(execFile);
const here = dirname(fileURLToPath(import.meta.url));
const packageRoot = dirname(here);

const ar = { model_id: "M1_AR", mode: "autosomal_recessive", target_label: "target" };
const assertFiniteTree = (value) => {
  if (typeof value === "number") assert.ok(Number.isFinite(value), "simulation must never emit non-finite numeric output");
  else if (Array.isArray(value)) value.forEach(assertFiniteTree);
  else if (value && typeof value === "object") Object.values(value).forEach(assertFiniteTree);
};
const hash = (character) => character.repeat(64);
const feedbackEnvelope = () => ({
  schema_version: "1.0.0", run_id: "CM-RUN-NODE-000001", operation: "cross", engine_version: "0.4.0", input_sha256: hash("a"), output_sha256: hash("b"),
  package: { name: "@biobirddna/genetics", version: "0.4.0", integrity_sha256: hash("c") },
  result: { evidence: { evidence_grade: "G0" }, limits: { decision_support_only: true, gene_status: "not_assessed" }, missing_data: [] }, human_review_status: "pending"
});

test("package metadata remains private-only", async () => {
  const packageJson = JSON.parse(await readFile(join(packageRoot, "package.json"), "utf8"));
  assert.equal(packageJson.private, true);
  for (const field of ["publishConfig", "repository", "bugs", "homepage"]) assert.equal(field in packageJson, false);
});

test("cross returns the deterministic AR distribution and limits", () => {
  const result = dispatch("cross", { model: ar, male_genotype: "N/m", female_genotype: "m/m" });
  assert.equal(result.phenotype_probabilities["male:target"], 0.25);
  assert.equal(result.phenotype_probabilities["female:target"], 0.25);
  assert.equal(result.limits.inheritance_status, "simulation_not_confirmation");
});

test("cross preserves the avian ZR reciprocal sex-specific distribution", () => {
  const result = dispatch("cross", {
    model: { model_id: "M4_ZR", mode: "z_linked_recessive", target_label: "target" },
    male_genotype: "m/m", female_genotype: "N/W"
  });
  assert.equal(result.phenotype_probabilities["male:normal"], 0.5);
  assert.equal(result.phenotype_probabilities["female:target"], 0.5);
});

test("cross rejects mixed input modes and non-genetic genotype aliases", () => {
  const m0 = { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" };
  assert.throws(
    () => dispatch("cross", { model: m0, male_genotype: "garbage", female_genotype: "NA" }),
    /non_genetic_model_requires_NA_genotype/
  );
  assert.throws(
    () => dispatch("cross", { model: m0, male_genotypes: { garbage: 1 }, female_genotypes: { NA: 1 } }),
    /non_genetic_model_requires_NA_genotype/
  );
  assert.throws(
    () => dispatch("cross", {
      model: ar,
      male_genotype: "N/m",
      female_genotype: "m/m",
      male_genotypes: { "N/m": 1 },
      female_genotypes: { "m/m": 1 }
    }),
    /cross_parent_input_modes_mutually_exclusive/
  );
  const valid = dispatch("cross", { model: m0, male_genotype: "NA", female_genotype: "NA" });
  assert.equal(valid.male_genotype, "NA");
  assert.equal(valid.female_genotype, "NA");
  assert.deepEqual(valid.limits, {
    gene_status: "not_assessed",
    interpretation: "non_genetic_baseline_model",
    claim_language: "simulate_cross"
  });
});

test("uncertain cross accepts nonnegative relative weights above one", () => {
  const result = dispatch("cross", {
    model: ar,
    male_genotypes: { "N/N": 2, "N/m": 1 },
    female_genotypes: { "m/m": 3 }
  });
  assert.ok(Math.abs(result.male_genotype_posterior["N/N"] - 2 / 3) < 1e-12);
  assert.deepEqual(result.female_genotype_posterior, { "m/m": 1 });
});

test("infer keeps posterior relative-only and includes M0", () => {
  const result = dispatch("infer", {
    target_label: "target",
    male_parent_phenotype: "target",
    female_parent_phenotype: "normal",
    offspring: [
      { sex: "male", phenotype: "normal", count: 5 },
      { sex: "female", phenotype: "target", count: 6 }
    ]
  });
  assert.equal(result.best_model, "M4_ZR");
  assert.ok(result.ranking.some((item) => item.model_id === "M0_NONGENETIC"));
  assert.equal(result.limits.gene_status, "not_assessed");
});

test("infer requires M0 plus another candidate model", () => {
  const request = {
    male_parent_phenotype: "normal",
    female_parent_phenotype: "normal",
    offspring: [{ sex: "unknown", phenotype: "normal", count: 1 }]
  };
  assert.throws(
    () => dispatch("infer", { ...request, models: [ar] }),
    /candidate_models_require_exactly_one_M0_NONGENETIC/
  );
  assert.throws(
    () => dispatch("infer", { ...request, models: [{ model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }] }),
    /candidate_models_require_multiple_models/
  );
  assert.throws(
    () => dispatch("infer", {
      ...request,
      models: [
        { model_id: "DUP", mode: "autosomal_recessive", target_label: "target" },
        { model_id: "DUP", mode: "autosomal_dominant", target_label: "target" },
        { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }
      ]
    }),
    /candidate_model_id_duplicate:DUP/
  );
});

test("direct infer output is JSON-safe and feedback-compatible", () => {
  const output = dispatch("infer", {
    target_label: "target",
    male_parent_phenotype: "unsupported",
    female_parent_phenotype: "unsupported",
    offspring: [{ sex: "unknown", phenotype: "normal", count: 1 }]
  });
  assertFiniteTree(output);
  assert.ok(output.ranking.some((item) => item.log_evidence === null));
  assert.equal(
    output.ranking
      .filter((item) => item.log_evidence === null)
      .reduce((sum, item) => sum + item.posterior, 0),
    0
  );
  assert.ok(Math.abs(output.ranking.reduce((sum, item) => sum + item.posterior, 0) - 1) < 1e-12);
  const envelope = feedbackEnvelope();
  envelope.operation = "infer";
  envelope.result = {
    output,
    evidence: { evidence_grade: "G0" },
    limits: { decision_support_only: true, gene_status: "not_assessed" },
    missing_data: []
  };
  assert.doesNotThrow(() => createFeedbackArtifact(envelope, hash("c")));
});

test("equal infer posteriors use deterministic model-id ordering", () => {
  const shared = { mode: "autosomal_recessive", target_label: "target" };
  const request = {
    male_parent_phenotype: "unknown",
    female_parent_phenotype: "unknown",
    offspring: [{ sex: "unknown", phenotype: "target", count: 2 }]
  };
  const m0 = { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" };
  const first = dispatch("infer", {
    ...request,
    models: [{ model_id: "Z_AR", ...shared }, { model_id: "A_AR", ...shared }, m0]
  });
  const second = dispatch("infer", {
    ...request,
    models: [{ model_id: "A_AR", ...shared }, { model_id: "Z_AR", ...shared }, m0]
  });
  assert.deepEqual(first, second);
  assert.equal(first.best_model, "A_AR");
  assert.deepEqual(first.ranking.slice(0, 2).map((item) => item.model_id), ["A_AR", "Z_AR"]);
});

test("infer defaults omitted offspring sex and rejects schema-invalid coercions", () => {
  const base = {
    target_label: "target",
    male_parent_phenotype: "normal",
    female_parent_phenotype: "normal",
    offspring: [{ phenotype: "normal", count: 1 }]
  };
  assert.doesNotThrow(() => dispatch("infer", base));
  for (const [change, error] of [
    [{ unknown: true }, /inference_request_unknown_field:unknown/],
    [{ target_label: 7 }, /target_label_required/],
    [{ male_parent_phenotype: 7 }, /male_parent_phenotype_required/],
    [{ offspring: [{ phenotype: "normal", count: "1" }] }, /offspring_count_invalid/],
    [{ offspring: [{ phenotype: "normal", count: true }] }, /offspring_count_invalid/],
    [{ offspring: [{ phenotype: "normal", count: 1, unknown: true }] }, /offspring_observation_unknown_field:unknown/]
  ]) assert.throws(() => dispatch("infer", { ...base, ...change }), error);
});

test("infer retains finite parent configurations when sibling configurations are impossible", () => {
  const result = dispatch("infer", {
    models: [
      { model_id: "M1_AR", mode: "autosomal_recessive", target_label: "target", observation_error: 0 },
      { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }
    ],
    male_parent_phenotype: "unknown",
    female_parent_phenotype: "unknown",
    offspring: [{ sex: "unknown", phenotype: "target", count: 1 }]
  });
  const arScore = result.ranking.find((item) => item.model_id === "M1_AR");
  assert.notEqual(arScore.log_evidence, null);
  assert.equal(arScore.configurations.length, 9);
  assert.ok(arScore.configurations.some((item) => item.log_likelihood === null));
  assert.ok(arScore.configurations.some((item) => Number.isFinite(item.log_likelihood)));
});

test("infer collapses unknown-sex offspring exactly as the Python runtime", () => {
  const result = dispatch("infer", { target_label: "target", male_parent_phenotype: "unknown", female_parent_phenotype: "unknown", offspring: [{ sex: "unknown", phenotype: "target", count: 10 }, { sex: "unknown", phenotype: "normal", count: 2 }] });
  const expected = { M5_ZD: [0.5275202877761334, -3.1214986246004037], M2_AD: [0.35793886528677676, -3.509323744987241], M4_ZR: [0.0674684028523687, -5.1780265670582795], M1_AR: [0.0449789352349125, -5.583491675166443], M3_INC: [0.0020935086866083015, -8.650844493658425], M0_NONGENETIC: [1.6320024985879716e-10, -25.01797381287021] };
  for (const item of result.ranking) { const [posterior, logEvidence] = expected[item.model_id]; assert.ok(Math.abs(item.posterior - posterior) < 2e-12, item.model_id); assert.ok(Math.abs(item.log_evidence - logEvidence) < 2e-12, item.model_id); }
});

test("pigment remains phenotype-only", () => {
  const result = dispatch("pigment", { white_area_ratio: 0.2, progressive_across_molts: true, image_quality: "Q3" });
  assert.equal(result.limits.gene_status, "not_assessed");
  assert.equal(result.limits.inheritance_status, "not_assessed");
  // Zero-score labels are not candidates: only MEL-PGR and MEL-LEU accrued score here.
  assert.deepEqual(result.alternative_hypotheses, ["MEL-LEU"]);
  for (const label of result.alternative_hypotheses) assert.ok(result.scores[label] > 0);
  // Control: all-positive ranking keeps the three-alternative shape.
  const allPositive = dispatch("pigment", { brown_shift_score: 0.6, grey_shift_score: 0.5, dark_expansion_score: 0.4, red_patch_change_score: 0.3, image_quality: "Q3" });
  assert.deepEqual(allPositive.alternative_hypotheses, ["MEL-DIL", "MEL-MEL", "CAR-ALT"]);
  for (const label of allPositive.alternative_hypotheses) assert.ok(allPositive.scores[label] > 0);
  for (const [request, error] of [
    [{ image_quality: "BOGUS" }, /image_quality_invalid/],
    [{ white_area_ratio: "0.5" }, /white_area_ratio_invalid/],
    [{ iris_red_confirmed: "false" }, /iris_red_confirmed_invalid/],
    [{ pattern_retained: "false" }, /pattern_retained_invalid/],
    [{ unknown: true }, /pigment_request_unknown_field:unknown/]
  ]) assert.throws(() => dispatch("pigment", request), error);
});

test("pigment forces Q0 and Q1 recapture despite strong visible features", () => {
  for (const image_quality of ["Q0", "Q1"]) {
    const result = dispatch("pigment", {
      white_area_ratio: 0.99,
      brown_shift_score: 1,
      grey_shift_score: 1,
      dark_expansion_score: 1,
      red_patch_change_score: 1,
      pattern_retained: true,
      iris_red_confirmed: true,
      progressive_across_molts: true,
      image_quality
    });
    assert.equal(result.primary_hypothesis, "UNK", image_quality);
    assert.equal(result.needs_recapture, true, image_quality);
    assert.equal(result.limits.measurement_domain, "visible_phenotype_only");
    assert.equal(result.limits.gene_status, "not_assessed");
    assert.equal(result.limits.inheritance_status, "not_assessed");
  }
});

test("rank-mates keeps a full-sibling pair blocked", () => {
  const result = dispatch("rank-mates", {
    pedigree: [
      { individual_id: "S", sex: "male", sire_id: null, dam_id: null },
      { individual_id: "D", sex: "female", sire_id: null, dam_id: null },
      { individual_id: "M", sex: "male", sire_id: "S", dam_id: "D" },
      { individual_id: "F", sex: "female", sire_id: "S", dam_id: "D" }
    ],
    males: [{ individual_id: "M", sex: "male", genotype_posteriors: { M1_AR: { "N/N": 1 } } }],
    females: [{ individual_id: "F", sex: "female", genotype_posteriors: { M1_AR: { "N/N": 1 } } }],
    hypotheses: [{ weight: 1, model: ar }]
  });
  assert.equal(result.pair_ranking[0].eligible, false);
  assert.ok(result.pair_ranking[0].block_reasons.includes("full_sibling_pair"));
});

test("rank-mates always blocks pairing an individual with itself", () => {
  const result = dispatch("rank-mates", {
    pedigree: [
      { individual_id: "SELF", sex: "male", sire_id: null, dam_id: null }
    ],
    males: [{ individual_id: "SELF", sex: "male", genotype_posteriors: { M1_AR: { "N/m": 1 } } }],
    females: [{ individual_id: "SELF", sex: "female", genotype_posteriors: { M1_AR: { "N/m": 1 } } }],
    hypotheses: [{ weight: 1, model: ar }],
    max_offspring_inbreeding: 1,
    block_direct_ancestor_and_full_sibling: false
  });
  assert.equal(result.pair_ranking[0].eligible, false);
  assert.ok(result.pair_ranking[0].block_reasons.includes("same_individual_pair"));
});

test("normalize-records preserves missing distinct from zero", () => {
  const result = dispatch("normalize-records", {
    schema_version: "1.0.0", case_id: "CM-CASE-NODE",
    pairings: [{ pairing_id: "P", male_parent_id: "M", female_parent_id: "F", male_parent_phenotype: "normal", female_parent_phenotype: "normal", male_parent_sex_method: "dna", male_parent_sex_confidence: 0.99, female_parent_sex_method: "dna", female_parent_sex_confidence: 0.99 }],
    clutches: [{ clutch_id: "C", pairing_id: "P", offspring_records_coverage: "complete", stage_counts: {
      eggs_laid: { status: "observed", coverage: "complete", count: 0 },
      fertile_eggs: { status: "observed", coverage: "complete", count: 0 },
      hatched: { status: "observed", coverage: "complete", count: 0 },
      survived_to_first_molt: { status: "missing", coverage: "unknown" }
    }}], offspring: []
  });
  assert.equal(result.clutches[0].stage_counts.eggs_laid.count, 0);
  assert.equal(result.clutches[0].stage_counts.survived_to_first_molt.count, null);
  assert.equal(result.evidence.genetics_grade, "not_assessed");
});

test("simulate-recovery is deterministic, G0, and requires M0", () => {
  const request = {
    schema_version: "1.0.0", base_seed: 7, goodness_of_fit_replicates: 20, misspecification_alpha: 0.05,
    experiments: [{ experiment_id: "SIM-AR", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [4], replicates: 2, candidate_model_ids: ["M1_AR", "M2_AD", "M3_INC", "M4_ZR", "M5_ZD", "M0_NONGENETIC"] }]
  };
  assert.deepEqual(dispatch("simulate-recovery", request), dispatch("simulate-recovery", request));
  assert.equal(dispatch("simulate-recovery", request).evidence.evidence_grade, "G0");
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], candidate_model_ids: ["M1_AR"] }] }), BioBirdDNAError);
});

test("simulate-recovery sampling does not require the Node Buffer global", () => {
  const bufferDescriptor = Object.getOwnPropertyDescriptor(globalThis, "Buffer");
  try {
    assert.equal(Reflect.deleteProperty(globalThis, "Buffer"), true);
    assert.equal(typeof globalThis.Buffer, "undefined");
    assert.deepEqual(__sampleForTest({ "𐀀": 1, "": 1 }, [0.25, 0.25]), { "": 2, "𐀀": 0 });
    assert.doesNotThrow(() => dispatch("simulate-recovery", {
      schema_version: "1.0.0", base_seed: 7, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1,
      experiments: [{ experiment_id: "SIM-BUFFER-FREE", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [2], replicates: 1, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }]
    }));
  } finally {
    if (bufferDescriptor) Object.defineProperty(globalThis, "Buffer", bufferDescriptor);
  }
});

test("simulate-recovery protects the GEN-011 preflight and bootstrap work gates", () => {
  const request = {
    schema_version: "1.0.0", base_seed: 9, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1,
    experiments: [{ experiment_id: "SIM-AR", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [4], replicates: 1, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }]
  };
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { ...ar, penetrance: 0.8 } }] }), /simulation_generating_model_snapshot_mismatch:SIM-AR:M1_AR/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { ...ar, model_id: null } }] }), /model_id_required/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { ...ar, model_id: " M1_AR " } }] }), /model_id_required/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { ...ar, target_label: null } }] }), /target_label_required/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { ...ar, penetrance: null } }] }), /penetrance_invalid/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], candidate_model_ids: [" M1_AR ", " M0_NONGENETIC "] }] }), /simulation_candidate_model_id_required/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { model_id: "M1_AR", mode: "autosomal_dominant", target_label: "target" }, candidate_model_ids: ["M2_AD", "M0_NONGENETIC"] }] }), /simulation_generating_model_snapshot_mismatch:SIM-AR:M1_AR/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], male_parent_phenotype: "target" }] }), /simulation_parent_phenotype_genotype_incompatible:SIM-AR:male/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }, male_genotype: "garbage", female_genotype: "NA", candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }] }), /simulation_non_genetic_genotype_must_be_NA:SIM-AR:male/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, experiments: [{ ...request.experiments[0], generating_model: { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }, male_genotype: "NA", female_genotype: "NA", candidate_model_ids: ["M0_NONGENETIC"] }] }), /simulation_candidate_set_requires_multiple_models:SIM-AR/);
  assert.throws(() => dispatch("simulate-recovery", { ...request, goodness_of_fit_replicates: 10000, experiments: [{ ...request.experiments[0], brood_sizes: [100], replicates: 100 }] }), /simulation_bootstrap_work_budget_exceeded/);
});

test("simulate-recovery rejects padded GEN-011 semantic strings and permits custom generators", () => {
  const base = { schema_version: "1.0.0", base_seed: 11, goodness_of_fit_replicates: 10, misspecification_alpha: .1, experiments: [{ experiment_id: "SIM-EXACT", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [2], replicates: 1, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }] };
  for (const change of [
    { experiment_id: " SIM-EXACT " }, { male_genotype: " N/m " }, { female_genotype: " N/m " }, { male_parent_phenotype: " normal " }, { female_parent_phenotype: " normal " }, { candidate_model_ids: [" M1_AR", "M0_NONGENETIC"] }, { generating_model: { ...ar, model_id: " M1_AR" } }, { generating_model: { ...ar, target_label: " target" } }
  ]) assert.throws(() => dispatch("simulate-recovery", { ...base, experiments: [{ ...base.experiments[0], ...change }] }), BioBirdDNAError);
  assert.doesNotThrow(() => dispatch("simulate-recovery", { ...base, experiments: [{ ...base.experiments[0], generating_model: { ...ar, model_id: "MISSPEC_AR", description: "custom misspecification" }, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }] }));
});

test("simulate-recovery preserves inference when phenotype labels contain colons", () => {
  const make = (target_label) => ({ schema_version: "1.0.0", base_seed: 21, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1, experiments: [{ experiment_id: "SIM-LABEL", generating_model: { ...ar, target_label }, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [8], replicates: 1, candidate_model_ids: ["M1_AR", "M2_AD", "M0_NONGENETIC"] }] });
  const plain = dispatch("simulate-recovery", make("target")).experiments[0].by_brood_size[0].replicate_results[0].model_posteriors;
  const colon = dispatch("simulate-recovery", make("target:shade")).experiments[0].by_brood_size[0].replicate_results[0].model_posteriors;
  assert.deepEqual(colon, plain);
});

test("simulate-recovery keeps every adversarial floating-weight draw in a bucket", () => {
  const result = dispatch("simulate-recovery", { schema_version: "1.0.0", base_seed: 123, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1, experiments: [{ experiment_id: "SIM-WEIGHTS", generating_model: { model_id: "SIM_NONGEN", mode: "non_genetic", target_label: "target", female_probability: 0.1, baseline_target_probability: 0.2, observation_error: 0.17 }, male_genotype: "NA", female_genotype: "NA", male_parent_phenotype: "unknown", female_parent_phenotype: "unknown", brood_sizes: [97], replicates: 3, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }] });
  const block = result.experiments[0].by_brood_size[0];
  for (const replicate of block.replicate_results) assert.equal(Object.values(replicate.observed_counts).reduce((sum, count) => sum + count, 0), 97);
  for (const candidate of block.goodness_of_fit.candidates) assert.ok(candidate.deviance === null || candidate.deviance >= 0);
});

test("simulate-recovery reports misspecification, tie-safe recovery, and safe deterministic seeds", () => {
  const base = {
    schema_version: "1.0.0", base_seed: 9007199254740991, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1,
    experiments: [{ experiment_id: "SIM-MISSPEC", generating_model: { model_id: "M3_INC", mode: "incomplete_dominant", target_label: "target" }, male_genotype: "N/N", female_genotype: "m/m", male_parent_phenotype: "normal", female_parent_phenotype: "target", brood_sizes: [40], replicates: 1, candidate_model_ids: ["M1_AR", "M2_AD", "M4_ZR", "M5_ZD", "M0_NONGENETIC"] }]
  };
  const misspec = dispatch("simulate-recovery", base).experiments[0].by_brood_size[0];
  assert.ok(misspec.warnings.includes("candidate_model_set_misspecified_suspected"));
  assert.ok(misspec.replicate_results[0].seed <= Number.MAX_SAFE_INTEGER);
  assertFiniteTree(misspec);
  const tieRequest = { ...base, experiments: [{ experiment_id: "SIM-TIE", generating_model: { model_id: "M2_AD", mode: "autosomal_dominant", target_label: "target" }, male_genotype: "N/N", female_genotype: "N/N", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [6], replicates: 2, candidate_model_ids: ["M2_AD", "M1_AR", "M3_INC", "M4_ZR", "M5_ZD", "M0_NONGENETIC"] }] };
  const reversed = { ...tieRequest, experiments: [{ ...tieRequest.experiments[0], candidate_model_ids: [...tieRequest.experiments[0].candidate_model_ids].reverse() }] };
  const first = dispatch("simulate-recovery", tieRequest), second = dispatch("simulate-recovery", reversed);
  assert.deepEqual(first.experiments[0].by_brood_size[0].recovery, second.experiments[0].by_brood_size[0].recovery);
  assert.ok(first.experiments[0].by_brood_size[0].replicate_results.some((result) => result.best_model === null && result.best_models.length > 1));
});

test("simulate-recovery calibrates only under its declared multi-model design", () => {
  const request = {
    schema_version: "1.0.0", base_seed: 11, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1,
    experiments: [
      { experiment_id: "SIM-AR", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [4], replicates: 1, candidate_model_ids: ["M1_AR", "M4_ZR", "M0_NONGENETIC"] },
      { experiment_id: "SIM-ZR", generating_model: { model_id: "M4_ZR", mode: "z_linked_recessive", target_label: "target" }, male_genotype: "m/m", female_genotype: "N/W", male_parent_phenotype: "target", female_parent_phenotype: "normal", brood_sizes: [4], replicates: 1, candidate_model_ids: ["M1_AR", "M4_ZR", "M0_NONGENETIC"] }
    ]
  };
  const result = dispatch("simulate-recovery", request);
  assert.equal(result.calibration_by_design[0].status, "estimated");
  assert.equal(result.calibration_by_design[0].method, "empirical_reliability_under_declared_simulation_design");
  assert.ok("observed_true_model_rate" in result.calibration_by_design[0].models[0].bins[0]);
});

test("feedback artifacts are deterministic, immutable, and provenance-bound", () => {
  const envelope = feedbackEnvelope(), first = createFeedbackArtifact(envelope, hash("c")), second = createFeedbackArtifact(structuredClone(envelope), hash("c"));
  assert.equal(first.schema_version, "biobirddna_feedback_v1");
  assert.equal(first.feedback_id, second.feedback_id);
  assert.equal(first.proposal.status, "pending_human_review");
  assert.equal(first.proposal.auto_apply, false);
  assert.deepEqual(first.proposal.target_files, []);
  assert.ok(Object.isFrozen(first));
  assert.throws(() => createFeedbackArtifact(envelope, hash("d")), /feedback_package_integrity_mismatch/);
});

test("feedback adapter fails closed for schema, mutation, and unsupported claim requests", () => {
  const envelope = feedbackEnvelope();
  assert.throws(() => createFeedbackArtifact({ ...envelope, unknown: true }, hash("c")), /feedback_envelope_unknown_field:unknown/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, mutation_request: "apply" } }, hash("c")), /feedback_mutation_request_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, engine_version: "0.2.1" }, hash("c")), /feedback_engine_version_mismatch/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, engine_version: " 0.4.0 " }, hash("c")), /feedback_engine_version_mismatch/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, input_sha256: ` ${hash("a")}` }, hash("c")), /feedback_input_sha256_invalid/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, limits: { decision_support_only: true, gene_status: "assessed" } } }, hash("c")), /feedback_gene_status_not_assessed_required/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, action: "apply" } }, hash("c")), /feedback_mutation_request_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, requested_action: "apply" } }, hash("c")), /feedback_mutation_request_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, proposal_status: "approved" } }, hash("c")), /feedback_mutation_request_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, output: { gene_status: "confirmed" } } }, hash("c")), /feedback_confirmed_gene_claim_requires_G5/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, output: { gene_status: "confirmed", finding: "causal gene identified" } } }, hash("c")), /feedback_confirmed_gene_claim_requires_G5/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, evidence: { evidence_grade: "G4", finding: "gene confirmed" } } }, hash("c")), /feedback_confirmed_gene_claim_requires_G5/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, evidence: { evidence_grade: "G4", finding: "confirmed gene" } } }, hash("c")), /feedback_confirmed_gene_claim_requires_G5/);
  for (const output of [
    { action: { type: "apply" } }, { apply: true }, { approval: true }, { write: { path: "src/index.mjs", content: "x" } }, { request: { kind: "delete" } },
    { gene: { status: "confirmed" } }, { claim: { subject: "gene", status: "identified" } }
  ]) assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, output } }, hash("c")));
});

test("feedback permits only explicit inert proposal markers", () => {
  const envelope = feedbackEnvelope();
  assert.doesNotThrow(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, target_files: [], auto_apply: false } }, hash("c")));
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, target_files: ["src/index.mjs"] } }, hash("c")), /feedback_mutation_request_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, auto_apply: true } }, hash("c")), /feedback_mutation_request_forbidden/);
});

test("feedback rejects numerical scoring fields at every nested level", () => {
  const envelope = feedbackEnvelope();
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, limits: { ...envelope.result.limits, health_gate_score: 3 } } }, hash("c")), /feedback_numerical_score_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, evidence: { evidence_grade: "G0", nested: { autoresearch_score: 0.2 } } } }, hash("c")), /feedback_numerical_score_forbidden/);
  assert.throws(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, output: { autoresearch: { health_gate: { score: 3.2 } } } } }, hash("c")), /feedback_numerical_score_forbidden/);
  assert.doesNotThrow(() => createFeedbackArtifact({ ...envelope, result: { ...envelope.result, fraction: 0.2, scores: { MEL: 0.8 } } }, hash("c")));
});

test("feedback accepts completed envelopes built from each public operation shape", () => {
  const normalize = dispatch("normalize-records", { schema_version: "1.0.0", case_id: "CM-FEEDBACK", pairings: [{ pairing_id: "P", male_parent_id: "M", female_parent_id: "F", male_parent_phenotype: "normal", female_parent_phenotype: "normal", male_parent_sex_method: "dna", male_parent_sex_confidence: 1, female_parent_sex_method: "dna", female_parent_sex_confidence: 1 }], clutches: [], offspring: [] });
  const simulation = dispatch("simulate-recovery", { schema_version: "1.0.0", base_seed: 5, goodness_of_fit_replicates: 10, misspecification_alpha: 0.1, experiments: [{ experiment_id: "SIM-FEEDBACK", generating_model: ar, male_genotype: "N/m", female_genotype: "N/m", male_parent_phenotype: "normal", female_parent_phenotype: "normal", brood_sizes: [2], replicates: 1, candidate_model_ids: ["M1_AR", "M0_NONGENETIC"] }] });
  const outputs = {
    cross: dispatch("cross", { model: ar, male_genotype: "N/m", female_genotype: "m/m" }),
    infer: dispatch("infer", { models: [ar, { model_id: "M0_NONGENETIC", mode: "non_genetic", target_label: "target" }], male_parent_phenotype: "normal", female_parent_phenotype: "normal", offspring: [{ sex: "male", phenotype: "normal", count: 1 }] }),
    pigment: dispatch("pigment", { image_quality: "Q3", white_area_ratio: 0.2 }),
    "rank-mates": dispatch("rank-mates", { pedigree: [{ individual_id: "M", sex: "male" }, { individual_id: "F", sex: "female" }], males: [{ individual_id: "M", sex: "male", genotype_posteriors: { M1_AR: { "N/N": 1 } } }], females: [{ individual_id: "F", sex: "female", genotype_posteriors: { M1_AR: { "N/N": 1 } } }], hypotheses: [{ model: ar, weight: 1 }] }),
    "normalize-records": normalize,
    "simulate-recovery": simulation
  };
  for (const [operation, output] of Object.entries(outputs)) {
    const envelope = { ...feedbackEnvelope(), operation, result: { output, evidence: output.evidence ?? {}, limits: { decision_support_only: true, gene_status: "not_assessed" }, missing_data: [] } };
    assert.doesNotThrow(() => createFeedbackArtifact(envelope, hash("c")), operation);
  }
});

test("feedback emits an insufficient-data signal without a numerical score", () => {
  const envelope = feedbackEnvelope(); envelope.result.evidence = {}; envelope.result.missing_data = ["phenotype observations absent"];
  const artifact = createFeedbackArtifact(envelope, hash("c"));
  assert.deepEqual(artifact.signals, [{ code: "INSUFFICIENT_DATA", severity: "warning", rationale: "Production evidence is missing or incomplete; no numerical score is emitted." }]);
  assert.equal(Object.keys(artifact).some((key) => /score/i.test(key)), false);
});

test("feedback CLI handles paths with spaces and adapter source has no prohibited imports", async () => {
  const dir = await mkdtemp(join(tmpdir(), "biobird feedback space "));
  try {
    const input = join(dir, "feedback input.json"), output = join(dir, "feedback output.json");
    await writeFile(input, JSON.stringify(feedbackEnvelope()));
    await execFileAsync(process.execPath, [join(packageRoot, "bin", "biobirddna-feedback.mjs"), "--input", input, "--expected-integrity", hash("c"), "--output", output]);
    assert.equal(JSON.parse(await readFile(output, "utf8")).schema_version, "biobirddna_feedback_v1");
    const source = await readFile(join(packageRoot, "src", "feedback.mjs"), "utf8");
    const cli = await readFile(join(packageRoot, "bin", "biobirddna-feedback.mjs"), "utf8");
    assert.doesNotMatch(source, /node:(?:http|https|net|tls|child_process)|AutoResearch|TAICHEL|fetch\s*\(/i);
    assert.ok(cli.indexOf("await stat") < cli.indexOf("await readFile"), "CLI must stat before reading input");
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("CLI works from a path containing spaces", async () => {
  const dir = await mkdtemp(join(tmpdir(), "biobird node space "));
  try {
    const input = join(dir, "input file.json");
    const output = join(dir, "output file.json");
    await writeFile(input, JSON.stringify({ model: ar, male_genotype: "N/m", female_genotype: "m/m" }));
    await execFileAsync(process.execPath, [join(packageRoot, "bin", "biobirddna-genetics.mjs"), "cross", "--input", input, "--output", output]);
    const result = JSON.parse(await readFile(output, "utf8"));
    assert.equal(result.ok, true);
    assert.equal(result.operation, "cross");
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});
