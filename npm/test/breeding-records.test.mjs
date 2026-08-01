import assert from "node:assert/strict";
import test from "node:test";

import {
  BreedingRecordValidationError,
  normalizeBreedingRecords
} from "../src/breeding-records.mjs";

function completeRecordRequest() {
  return {
    schema_version: "1.0.0",
    case_id: "CM-CASE-000001",
    pairings: [
      {
        pairing_id: "CM-MAT-000002",
        male_parent_id: "CM-IND-B",
        female_parent_id: "CM-IND-A",
        male_parent_phenotype: "target_colour",
        female_parent_phenotype: "normal",
        male_parent_sex_method: "dna",
        male_parent_sex_confidence: 0.99,
        female_parent_sex_method: "dna",
        female_parent_sex_confidence: 0.99,
        reciprocal_group_id: "CM-RECIP-000001"
      },
      {
        pairing_id: "CM-MAT-000001",
        male_parent_id: "CM-IND-A",
        female_parent_id: "CM-IND-B",
        male_parent_phenotype: "normal",
        female_parent_phenotype: "target_colour",
        male_parent_sex_method: "dna",
        male_parent_sex_confidence: 0.99,
        female_parent_sex_method: "dna",
        female_parent_sex_confidence: 0.99,
        reciprocal_group_id: "CM-RECIP-000001"
      }
    ],
    clutches: [
      {
        clutch_id: "CM-CLU-000001",
        pairing_id: "CM-MAT-000001",
        offspring_records_coverage: "complete",
        stage_counts: {
          eggs_laid: { status: "observed", coverage: "complete", count: 4 },
          fertile_eggs: { status: "observed", coverage: "complete", count: 3 },
          hatched: { status: "observed", coverage: "complete", count: 2 },
          survived_to_first_molt: { status: "missing", coverage: "unknown" }
        }
      },
      {
        clutch_id: "CM-CLU-000002",
        pairing_id: "CM-MAT-000002",
        offspring_records_coverage: "complete",
        stage_counts: {
          eggs_laid: { status: "observed", coverage: "complete", count: 0 },
          fertile_eggs: { status: "observed", coverage: "complete", count: 0 },
          hatched: { status: "observed", coverage: "complete", count: 0 },
          survived_to_first_molt: {
            status: "observed",
            coverage: "complete",
            count: 0
          }
        }
      }
    ],
    offspring: [
      {
        offspring_record_id: "CM-OFF-000001",
        individual_id: "CM-IND-OFF-000001",
        clutch_id: "CM-CLU-000001",
        outcome: "alive",
        fertility_status: "fertile",
        hatch_status: "hatched",
        survival_status: "survived",
        sex: "male",
        sex_method: "dna",
        sex_confidence: 0.99,
        phenotype_timeline: [
          {
            observation_id: "CM-PHENO-000001",
            sequence_index: 1,
            stage: "first_plumage",
            status: "observed",
            phenotype: "normal"
          }
        ]
      },
      {
        offspring_record_id: "CM-OFF-000002",
        individual_id: "CM-IND-OFF-000002",
        clutch_id: "CM-CLU-000001",
        outcome: "post_hatch_death",
        fertility_status: "fertile",
        hatch_status: "hatched",
        survival_status: "died_before_first_molt",
        sex: "unknown",
        sex_method: "unknown",
        sex_confidence: null,
        phenotype_timeline: [
          {
            observation_id: "CM-PHENO-000002",
            sequence_index: 0,
            stage: "hatch",
            status: "missing"
          }
        ]
      },
      {
        offspring_record_id: "CM-OFF-000003",
        individual_id: null,
        clutch_id: "CM-CLU-000001",
        outcome: "embryo_death",
        fertility_status: "fertile",
        hatch_status: "not_hatched",
        survival_status: "not_applicable",
        sex: "unknown",
        sex_method: "unknown",
        sex_confidence: null,
        phenotype_timeline: []
      },
      {
        offspring_record_id: "CM-OFF-000004",
        individual_id: null,
        clutch_id: "CM-CLU-000001",
        outcome: "unfertilized",
        fertility_status: "unfertilized",
        hatch_status: "not_applicable",
        survival_status: "not_applicable",
        sex: "unknown",
        sex_method: "unknown",
        sex_confidence: null,
        phenotype_timeline: []
      }
    ]
  };
}

function assertValidation(request, expectedMessage) {
  assert.throws(
    () => normalizeBreedingRecords(request),
    (error) => {
      assert.ok(error instanceof BreedingRecordValidationError);
      assert.equal(error.name, "BreedingRecordValidationError");
      assert.equal(error.code, expectedMessage);
      assert.equal(error.message, expectedMessage);
      return true;
    }
  );
}

test("normalizes GEN-010 records without mutation and preserves reciprocal direction", () => {
  const request = completeRecordRequest();
  const original = structuredClone(request);

  const result = normalizeBreedingRecords(request);

  assert.deepEqual(request, original);
  assert.deepEqual(
    result.pairings.map(({ pairing_id, cross_direction_key }) => ({
      pairing_id,
      cross_direction_key
    })),
    [
      {
        pairing_id: "CM-MAT-000001",
        cross_direction_key: "CM-IND-A>CM-IND-B"
      },
      {
        pairing_id: "CM-MAT-000002",
        cross_direction_key: "CM-IND-B>CM-IND-A"
      }
    ]
  );
  assert.deepEqual(
    result.offspring.map(({ outcome }) => outcome),
    ["alive", "post_hatch_death", "embryo_death", "unfertilized"]
  );
  assert.equal(
    result.offspring[0].phenotype_timeline[0].phenotype,
    "normal"
  );
  assert.deepEqual(result.evidence, {
    source: "submitted_breeding_records",
    genetics_grade: "not_assessed",
    record_scope: "normalized_observations_only"
  });
  assert.deepEqual(result.limits, {
    purpose: "breeding_record_normalization_only",
    gene_status: "not_assessed",
    inheritance_status: "not_assessed",
    decision_support_only: true,
    sex_method_recorded_not_verified: true
  });
});

test("keeps missing stages distinct from observed zero and reports missing data", () => {
  const result = normalizeBreedingRecords(completeRecordRequest());

  assert.deepEqual(
    result.clutches[0].stage_counts.survived_to_first_molt,
    { status: "missing", coverage: "unknown", count: null }
  );
  assert.deepEqual(result.clutches[1].stage_counts.eggs_laid, {
    status: "observed",
    coverage: "complete",
    count: 0
  });
  assert.deepEqual(result.missing_data, [
    {
      field:
        "clutches/CM-CLU-000001/stage_counts/survived_to_first_molt",
      impact: "survival_stage_not_observed"
    },
    {
      field:
        "offspring/CM-OFF-000002/phenotype_timeline/CM-PHENO-000002",
      impact: "phenotype_stage_not_observed"
    }
  ]);
});

test("sorts phenotype observations and rejects duplicate timeline identity", () => {
  const request = completeRecordRequest();
  request.offspring[0].phenotype_timeline.push(
    {
      observation_id: "CM-PHENO-000003",
      sequence_index: 2,
      stage: "adult",
      status: "observed",
      phenotype: "normal"
    },
    {
      observation_id: "CM-PHENO-000004",
      sequence_index: 0,
      stage: "hatch",
      status: "missing"
    }
  );
  assert.deepEqual(
    normalizeBreedingRecords(request).offspring[0].phenotype_timeline.map(
      ({ sequence_index }) => sequence_index
    ),
    [0, 1, 2]
  );

  const duplicateSequence = completeRecordRequest();
  duplicateSequence.offspring[0].phenotype_timeline.push({
    observation_id: "CM-PHENO-000003",
    sequence_index: 1,
    stage: "adult",
    status: "observed",
    phenotype: "normal"
  });
  assertValidation(
    duplicateSequence,
    "phenotype_sequence_index_duplicate:CM-OFF-000001:1"
  );

  const duplicateId = completeRecordRequest();
  duplicateId.offspring[1].phenotype_timeline[0].observation_id =
    "CM-PHENO-000001";
  assertValidation(
    duplicateId,
    "phenotype_observation_id_duplicate:CM-PHENO-000001"
  );
});

test("rejects unknown fields at every object level", () => {
  const cases = [
    [
      (request) => {
        request.inferred_gene = "forbidden";
      },
      "breeding_record_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.pairings[0].inferred_gene = "forbidden";
      },
      "pairing_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.clutches[0].inferred_gene = "forbidden";
      },
      "clutch_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.clutches[0].stage_counts.inferred_gene = {
          status: "missing",
          coverage: "unknown"
        };
      },
      "stage_counts_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.clutches[0].stage_counts.eggs_laid.inferred_gene =
          "forbidden";
      },
      "stage_count_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.offspring[0].inferred_gene = "forbidden";
      },
      "offspring_record_unknown_field:inferred_gene"
    ],
    [
      (request) => {
        request.offspring[0].phenotype_timeline[0].inferred_gene =
          "forbidden";
      },
      "phenotype_observation_unknown_field:inferred_gene"
    ]
  ];

  for (const [mutate, expected] of cases) {
    const request = completeRecordRequest();
    mutate(request);
    assertValidation(request, expected);
  }
});

test("requires nonempty pairings and explicit parent sex provenance", () => {
  const empty = completeRecordRequest();
  empty.pairings = [];
  assertValidation(empty, "pairings_must_not_be_empty");

  const unknownMethod = completeRecordRequest();
  unknownMethod.pairings[0].male_parent_sex_method = "unknown";
  unknownMethod.pairings[0].male_parent_sex_confidence = null;
  assertValidation(
    unknownMethod,
    "parent_sex_provenance_required:CM-MAT-000002:male"
  );

  const invalidConfidence = completeRecordRequest();
  invalidConfidence.pairings[0].female_parent_sex_confidence = Number.NaN;
  assertValidation(
    invalidConfidence,
    "parent_sex_provenance_required:CM-MAT-000002:female"
  );
});

test("requires every explicit stage and validates missing versus observed status", () => {
  const absentStage = completeRecordRequest();
  delete absentStage.clutches[0].stage_counts.hatched;
  assertValidation(
    absentStage,
    "stage_count_required:CM-CLU-000001:hatched"
  );

  const missingWithCount = completeRecordRequest();
  missingWithCount.clutches[0].stage_counts.survived_to_first_molt.count = 0;
  assertValidation(
    missingWithCount,
    "missing_stage_count_must_be_null:CM-CLU-000001:survived_to_first_molt"
  );

  const observedWithoutInteger = completeRecordRequest();
  observedWithoutInteger.clutches[0].stage_counts.eggs_laid.count = 4.5;
  assertValidation(
    observedWithoutInteger,
    "stage_count_invalid:CM-CLU-000001:eggs_laid"
  );
});

test("records partial coverage and reconciles only complete observations", () => {
  const partialStage = completeRecordRequest();
  partialStage.clutches[0].stage_counts.hatched.coverage = "partial";
  partialStage.clutches[0].stage_counts.hatched.count = 99;
  const result = normalizeBreedingRecords(partialStage);
  assert.ok(
    result.missing_data.some(
      (item) =>
        item.field === "clutches/CM-CLU-000001/stage_counts/hatched" &&
        item.impact === "hatch_stage_incomplete_coverage"
    )
  );

  const conflict = completeRecordRequest();
  conflict.clutches[0].stage_counts.eggs_laid.count = 3;
  assertValidation(
    conflict,
    "stage_offspring_count_conflict:CM-CLU-000001:eggs_laid:3:4"
  );

  const increasing = completeRecordRequest();
  increasing.clutches[0].stage_counts.hatched.count = 4;
  assertValidation(
    increasing,
    "clutch_stage_count_increases:CM-CLU-000001:fertile_eggs:3:hatched:4"
  );
});

test("validates outcome paths and refuses unknown status in complete reconciliation", () => {
  const conflict = completeRecordRequest();
  Object.assign(conflict.offspring[0], {
    fertility_status: "unfertilized",
    hatch_status: "not_applicable",
    survival_status: "not_applicable"
  });
  assertValidation(
    conflict,
    "offspring_outcome_stage_conflict:CM-OFF-000001:alive"
  );

  const unknownStatus = completeRecordRequest();
  Object.assign(unknownStatus.offspring[0], {
    outcome: "lost_to_followup",
    fertility_status: "fertile",
    hatch_status: "unknown",
    survival_status: "unknown"
  });
  assertValidation(
    unknownStatus,
    "offspring_stage_status_unknown:CM-CLU-000001:hatched:CM-OFF-000001"
  );
});

test("requires offspring provenance fields and method for a known sex", () => {
  for (const field of ["sex_confidence", "phenotype_timeline"]) {
    const request = completeRecordRequest();
    delete request.offspring[1][field];
    assertValidation(
      request,
      `offspring_field_required:CM-OFF-000002:${field}`
    );
  }

  const knownSex = completeRecordRequest();
  knownSex.offspring[0].sex_method = "unknown";
  knownSex.offspring[0].sex_confidence = null;
  assertValidation(
    knownSex,
    "known_sex_requires_method_and_confidence:CM-OFF-000001"
  );
});

test("validates schema identity, references, and globally unique record IDs", () => {
  const schema = completeRecordRequest();
  schema.schema_version = "2.0.0";
  assertValidation(schema, "breeding_schema_version_unsupported:2.0.0");

  const unknownPairing = completeRecordRequest();
  unknownPairing.clutches[0].pairing_id = "missing";
  assertValidation(
    unknownPairing,
    "clutch_pairing_unknown:CM-CLU-000001:missing"
  );

  const unknownClutch = completeRecordRequest();
  unknownClutch.offspring[0].clutch_id = "missing";
  assertValidation(
    unknownClutch,
    "offspring_clutch_unknown:CM-OFF-000001:missing"
  );

  const duplicate = completeRecordRequest();
  duplicate.offspring[1].offspring_record_id = "CM-OFF-000001";
  assertValidation(
    duplicate,
    "offspring_record_id_duplicate:CM-OFF-000001"
  );
});
