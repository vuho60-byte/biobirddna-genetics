const STAGE_ORDER = [
  "eggs_laid",
  "fertile_eggs",
  "hatched",
  "survived_to_first_molt"
];

const STAGE_MISSING_IMPACT = {
  eggs_laid: "egg_stage_not_observed",
  fertile_eggs: "fertility_stage_not_observed",
  hatched: "hatch_stage_not_observed",
  survived_to_first_molt: "survival_stage_not_observed"
};

const STAGE_SCOPE = {
  eggs_laid: "egg",
  fertile_eggs: "fertility",
  hatched: "hatch",
  survived_to_first_molt: "survival"
};

const OUTCOMES = new Set([
  "alive",
  "unfertilized",
  "embryo_death",
  "unhatched",
  "hatch_failure",
  "post_hatch_death",
  "lost_to_followup",
  "unknown"
]);
const SEXES = new Set(["male", "female", "unknown"]);
const SEX_METHODS = new Set([
  "dna",
  "reproductive",
  "veterinary",
  "visual",
  "reported",
  "unknown"
]);
const VERIFIED_SEX_METHODS = new Set(
  [...SEX_METHODS].filter((method) => method !== "unknown")
);
const FERTILITY_STATUSES = new Set(["fertile", "unfertilized", "unknown"]);
const HATCH_STATUSES = new Set([
  "hatched",
  "not_hatched",
  "not_applicable",
  "unknown"
]);
const SURVIVAL_STATUSES = new Set([
  "survived",
  "died_before_first_molt",
  "not_applicable",
  "unknown"
]);

function pathKey(fertilityStatus, hatchStatus, survivalStatus) {
  return `${fertilityStatus}\u0000${hatchStatus}\u0000${survivalStatus}`;
}

const VALID_STAGE_PATHS = new Set([
  pathKey("unfertilized", "not_applicable", "not_applicable"),
  pathKey("fertile", "not_hatched", "not_applicable"),
  pathKey("fertile", "hatched", "survived"),
  pathKey("fertile", "hatched", "died_before_first_molt"),
  pathKey("fertile", "hatched", "unknown"),
  pathKey("fertile", "unknown", "unknown"),
  pathKey("unknown", "unknown", "unknown")
]);

const OUTCOME_PATHS = {
  alive: new Set([
    pathKey("fertile", "hatched", "survived"),
    pathKey("fertile", "hatched", "unknown")
  ]),
  unfertilized: new Set([
    pathKey("unfertilized", "not_applicable", "not_applicable")
  ]),
  embryo_death: new Set([
    pathKey("fertile", "not_hatched", "not_applicable")
  ]),
  unhatched: new Set([
    pathKey("fertile", "not_hatched", "not_applicable")
  ]),
  hatch_failure: new Set([
    pathKey("fertile", "not_hatched", "not_applicable")
  ]),
  post_hatch_death: new Set([
    pathKey("fertile", "hatched", "died_before_first_molt")
  ]),
  lost_to_followup: new Set([
    pathKey("fertile", "hatched", "unknown"),
    pathKey("fertile", "unknown", "unknown"),
    pathKey("unknown", "unknown", "unknown")
  ]),
  unknown: VALID_STAGE_PATHS
};

const TOP_LEVEL_FIELDS = new Set([
  "schema_version",
  "case_id",
  "pairings",
  "clutches",
  "offspring"
]);
const PAIRING_FIELDS = new Set([
  "pairing_id",
  "male_parent_id",
  "female_parent_id",
  "male_parent_phenotype",
  "female_parent_phenotype",
  "male_parent_sex_method",
  "male_parent_sex_confidence",
  "female_parent_sex_method",
  "female_parent_sex_confidence",
  "reciprocal_group_id"
]);
const CLUTCH_FIELDS = new Set([
  "clutch_id",
  "pairing_id",
  "offspring_records_coverage",
  "stage_counts"
]);
const STAGE_COUNT_FIELDS = new Set(["status", "coverage", "count"]);
const PHENOTYPE_OBSERVATION_FIELDS = new Set([
  "observation_id",
  "sequence_index",
  "stage",
  "status",
  "phenotype"
]);
const OFFSPRING_FIELDS = new Set([
  "offspring_record_id",
  "individual_id",
  "clutch_id",
  "outcome",
  "fertility_status",
  "hatch_status",
  "survival_status",
  "sex",
  "sex_method",
  "sex_confidence",
  "phenotype_timeline"
]);

export class BreedingRecordValidationError extends Error {
  constructor(code) {
    super(code);
    this.name = "BreedingRecordValidationError";
    this.code = code;
  }
}

export function adaptBreedingRecordError(error, ErrorConstructor = Error) {
  if (!(error instanceof BreedingRecordValidationError)) {
    return error;
  }
  return new ErrorConstructor(error.message);
}

function fail(code) {
  throw new BreedingRecordValidationError(code);
}

function mapping(value, name) {
  if (value === null || typeof value !== "object" || Array.isArray(value)) {
    fail(`${name}_must_be_object`);
  }
  return value;
}

function sequence(value, name) {
  if (!Array.isArray(value)) {
    fail(`${name}_must_be_array`);
  }
  return value;
}

function requiredText(value, name) {
  if (typeof value !== "string" || value.trim().length === 0) {
    fail(`${name}_required`);
  }
  return value;
}

function optionalText(value, name) {
  if (value === null || value === undefined) {
    return null;
  }
  return requiredText(value, name);
}

function rejectUnknownFields(item, allowed, name) {
  const unknown = Object.keys(item)
    .filter((field) => !allowed.has(field))
    .sort();
  if (unknown.length > 0) {
    fail(`${name}_unknown_field:${unknown[0]}`);
  }
}

function uniqueId(value, name, seen) {
  const identifier = requiredText(value, name);
  if (seen.has(identifier)) {
    fail(`${name}_duplicate:${identifier}`);
  }
  seen.add(identifier);
  return identifier;
}

function normalizeStageCount(value, clutchId, stage, missingData) {
  const item = mapping(value, `${stage}_stage_count`);
  rejectUnknownFields(item, STAGE_COUNT_FIELDS, "stage_count");
  const status = item.status;
  if (status === "missing") {
    if (item.coverage !== "unknown") {
      fail(`missing_stage_coverage_must_be_unknown:${clutchId}:${stage}`);
    }
    if (item.count !== null && item.count !== undefined) {
      fail(`missing_stage_count_must_be_null:${clutchId}:${stage}`);
    }
    missingData.push({
      field: `clutches/${clutchId}/stage_counts/${stage}`,
      impact: STAGE_MISSING_IMPACT[stage]
    });
    return { status: "missing", coverage: "unknown", count: null };
  }
  if (status !== "observed") {
    fail(`stage_count_status_invalid:${clutchId}:${stage}`);
  }
  const coverage = item.coverage;
  if (coverage !== "complete" && coverage !== "partial") {
    fail(`stage_coverage_invalid:${clutchId}:${stage}`);
  }
  const count = item.count;
  if (!Number.isInteger(count) || count < 0) {
    fail(`stage_count_invalid:${clutchId}:${stage}`);
  }
  if (coverage === "partial") {
    missingData.push({
      field: `clutches/${clutchId}/stage_counts/${stage}`,
      impact: `${STAGE_SCOPE[stage]}_stage_incomplete_coverage`
    });
  }
  return { status: "observed", coverage, count };
}

function normalizePairings(value) {
  const pairings = [];
  const pairingIds = new Set();
  const rawPairings = sequence(value, "pairings");
  if (rawPairings.length === 0) {
    fail("pairings_must_not_be_empty");
  }
  for (const raw of rawPairings) {
    const item = mapping(raw, "pairing");
    rejectUnknownFields(item, PAIRING_FIELDS, "pairing");
    const pairingId = uniqueId(item.pairing_id, "pairing_id", pairingIds);
    const maleParentId = requiredText(item.male_parent_id, "male_parent_id");
    const femaleParentId = requiredText(
      item.female_parent_id,
      "female_parent_id"
    );
    if (maleParentId === femaleParentId) {
      fail(`pairing_parents_must_differ:${pairingId}`);
    }
    const parentSexProvenance = {};
    for (const role of ["male", "female"]) {
      const method = item[`${role}_parent_sex_method`];
      const confidence = item[`${role}_parent_sex_confidence`];
      if (
        !VERIFIED_SEX_METHODS.has(method) ||
        typeof confidence !== "number" ||
        !Number.isFinite(confidence) ||
        confidence < 0 ||
        confidence > 1
      ) {
        fail(`parent_sex_provenance_required:${pairingId}:${role}`);
      }
      parentSexProvenance[`${role}_parent_sex_method`] = method;
      parentSexProvenance[`${role}_parent_sex_confidence`] = confidence;
    }
    pairings.push({
      pairing_id: pairingId,
      male_parent_id: maleParentId,
      female_parent_id: femaleParentId,
      male_parent_phenotype: requiredText(
        item.male_parent_phenotype,
        "male_parent_phenotype"
      ),
      female_parent_phenotype: requiredText(
        item.female_parent_phenotype,
        "female_parent_phenotype"
      ),
      ...parentSexProvenance,
      reciprocal_group_id: optionalText(
        item.reciprocal_group_id,
        "reciprocal_group_id"
      ),
      cross_direction_key: `${maleParentId}>${femaleParentId}`
    });
  }
  pairings.sort((left, right) =>
    left.pairing_id < right.pairing_id
      ? -1
      : left.pairing_id > right.pairing_id
        ? 1
        : 0
  );
  return { pairings, pairingIds };
}

function normalizeClutches(value, pairingIds, missingData) {
  const clutches = [];
  const clutchIds = new Set();
  for (const raw of sequence(value, "clutches")) {
    const item = mapping(raw, "clutch");
    rejectUnknownFields(item, CLUTCH_FIELDS, "clutch");
    const clutchId = uniqueId(item.clutch_id, "clutch_id", clutchIds);
    const pairingId = requiredText(item.pairing_id, "pairing_id");
    if (!pairingIds.has(pairingId)) {
      fail(`clutch_pairing_unknown:${clutchId}:${pairingId}`);
    }
    const recordsCoverage = item.offspring_records_coverage;
    if (
      recordsCoverage !== "complete" &&
      recordsCoverage !== "partial" &&
      recordsCoverage !== "unknown"
    ) {
      fail(`offspring_records_coverage_invalid:${clutchId}`);
    }
    if (recordsCoverage !== "complete") {
      missingData.push({
        field: `clutches/${clutchId}/offspring`,
        impact:
          recordsCoverage === "partial"
            ? "offspring_records_incomplete_coverage"
            : "offspring_records_unknown_coverage"
      });
    }
    const rawCounts = mapping(item.stage_counts, "stage_counts");
    rejectUnknownFields(rawCounts, new Set(STAGE_ORDER), "stage_counts");
    const stageCounts = {};
    for (const stage of STAGE_ORDER) {
      if (!Object.hasOwn(rawCounts, stage)) {
        fail(`stage_count_required:${clutchId}:${stage}`);
      }
      stageCounts[stage] = normalizeStageCount(
        rawCounts[stage],
        clutchId,
        stage,
        missingData
      );
    }
    let previousStage = null;
    let previousCount = null;
    for (const stage of STAGE_ORDER) {
      const observation = stageCounts[stage];
      if (
        observation.status !== "observed" ||
        observation.coverage !== "complete"
      ) {
        continue;
      }
      const count = observation.count;
      if (previousCount !== null && count > previousCount) {
        fail(
          `clutch_stage_count_increases:${clutchId}:${previousStage}:` +
            `${previousCount}:${stage}:${count}`
        );
      }
      previousStage = stage;
      previousCount = count;
    }
    clutches.push({
      clutch_id: clutchId,
      pairing_id: pairingId,
      offspring_records_coverage: recordsCoverage,
      stage_counts: stageCounts
    });
  }
  clutches.sort((left, right) =>
    left.clutch_id < right.clutch_id
      ? -1
      : left.clutch_id > right.clutch_id
        ? 1
        : 0
  );
  return { clutches, clutchIds };
}

function confidence(value, offspringRecordId) {
  if (value === null || value === undefined) {
    return null;
  }
  if (
    typeof value !== "number" ||
    !Number.isFinite(value) ||
    value < 0 ||
    value > 1
  ) {
    fail(`sex_confidence_invalid:${offspringRecordId}`);
  }
  return value;
}

function normalizeTimeline(
  value,
  offspringRecordId,
  seenObservationIds,
  missingData
) {
  const timeline = [];
  const sequenceIndices = new Set();
  for (const raw of sequence(value, "phenotype_timeline")) {
    const item = mapping(raw, "phenotype_observation");
    rejectUnknownFields(
      item,
      PHENOTYPE_OBSERVATION_FIELDS,
      "phenotype_observation"
    );
    const observationId = uniqueId(
      item.observation_id,
      "phenotype_observation_id",
      seenObservationIds
    );
    const sequenceIndex = item.sequence_index;
    if (!Number.isInteger(sequenceIndex) || sequenceIndex < 0) {
      fail(`phenotype_sequence_index_invalid:${observationId}`);
    }
    if (sequenceIndices.has(sequenceIndex)) {
      fail(
        `phenotype_sequence_index_duplicate:${offspringRecordId}:` +
          `${sequenceIndex}`
      );
    }
    sequenceIndices.add(sequenceIndex);
    const stage = requiredText(item.stage, "phenotype_stage");
    const status = item.status;
    let phenotype;
    if (status === "observed") {
      phenotype = requiredText(item.phenotype, "phenotype");
    } else if (status === "missing") {
      if (item.phenotype !== null && item.phenotype !== undefined) {
        fail(`missing_phenotype_must_be_null:${observationId}`);
      }
      phenotype = null;
      missingData.push({
        field:
          `offspring/${offspringRecordId}/phenotype_timeline/` +
          observationId,
        impact: "phenotype_stage_not_observed"
      });
    } else {
      fail(`phenotype_status_invalid:${observationId}`);
    }
    timeline.push({
      observation_id: observationId,
      sequence_index: sequenceIndex,
      stage,
      status,
      phenotype
    });
  }
  timeline.sort(
    (left, right) => left.sequence_index - right.sequence_index
  );
  return timeline;
}

function stageStatus(value, allowed, name, offspringRecordId) {
  if (!allowed.has(value)) {
    fail(`${name}_invalid:${offspringRecordId}`);
  }
  return value;
}

function validateStagePath(
  offspringRecordId,
  outcome,
  fertilityStatus,
  hatchStatus,
  survivalStatus
) {
  const stagePath = pathKey(fertilityStatus, hatchStatus, survivalStatus);
  if (!VALID_STAGE_PATHS.has(stagePath)) {
    fail(`offspring_stage_path_invalid:${offspringRecordId}`);
  }
  if (!OUTCOME_PATHS[outcome].has(stagePath)) {
    fail(`offspring_outcome_stage_conflict:${offspringRecordId}:${outcome}`);
  }
}

function normalizeOffspring(value, clutchIds, missingData) {
  const offspring = [];
  const offspringIds = new Set();
  const observationIds = new Set();
  for (const raw of sequence(value, "offspring")) {
    const item = mapping(raw, "offspring_record");
    rejectUnknownFields(item, OFFSPRING_FIELDS, "offspring_record");
    const offspringRecordId = uniqueId(
      item.offspring_record_id,
      "offspring_record_id",
      offspringIds
    );
    const clutchId = requiredText(item.clutch_id, "clutch_id");
    if (!clutchIds.has(clutchId)) {
      fail(`offspring_clutch_unknown:${offspringRecordId}:${clutchId}`);
    }
    const outcome = item.outcome;
    if (!OUTCOMES.has(outcome)) {
      fail(`offspring_outcome_invalid:${offspringRecordId}`);
    }
    for (const requiredField of ["sex_confidence", "phenotype_timeline"]) {
      if (!Object.hasOwn(item, requiredField)) {
        fail(`offspring_field_required:${offspringRecordId}:${requiredField}`);
      }
    }
    const fertilityStatus = stageStatus(
      item.fertility_status,
      FERTILITY_STATUSES,
      "fertility_status",
      offspringRecordId
    );
    const hatchStatus = stageStatus(
      item.hatch_status,
      HATCH_STATUSES,
      "hatch_status",
      offspringRecordId
    );
    const survivalStatus = stageStatus(
      item.survival_status,
      SURVIVAL_STATUSES,
      "survival_status",
      offspringRecordId
    );
    validateStagePath(
      offspringRecordId,
      outcome,
      fertilityStatus,
      hatchStatus,
      survivalStatus
    );
    const sex = item.sex;
    if (!SEXES.has(sex)) {
      fail(`offspring_sex_invalid:${offspringRecordId}`);
    }
    const sexMethod = item.sex_method;
    if (!SEX_METHODS.has(sexMethod)) {
      fail(`sex_method_invalid:${offspringRecordId}`);
    }
    const sexConfidence = confidence(item.sex_confidence, offspringRecordId);
    if (
      sex !== "unknown" &&
      (sexMethod === "unknown" || sexConfidence === null)
    ) {
      fail(`known_sex_requires_method_and_confidence:${offspringRecordId}`);
    }
    offspring.push({
      offspring_record_id: offspringRecordId,
      individual_id: optionalText(item.individual_id, "individual_id"),
      clutch_id: clutchId,
      outcome,
      fertility_status: fertilityStatus,
      hatch_status: hatchStatus,
      survival_status: survivalStatus,
      sex,
      sex_method: sexMethod,
      sex_confidence: sexConfidence,
      phenotype_timeline: normalizeTimeline(
        item.phenotype_timeline,
        offspringRecordId,
        observationIds,
        missingData
      )
    });
  }
  return offspring;
}

function reconcileCompleteCounts(clutches, offspring) {
  const recordsByClutch = Object.fromEntries(
    clutches.map((clutch) => [clutch.clutch_id, []])
  );
  for (const record of offspring) {
    recordsByClutch[record.clutch_id].push(record);
  }
  const statusFields = {
    fertile_eggs: ["fertility_status", "fertile"],
    hatched: ["hatch_status", "hatched"],
    survived_to_first_molt: ["survival_status", "survived"]
  };
  for (const clutch of clutches) {
    if (clutch.offspring_records_coverage !== "complete") {
      continue;
    }
    const clutchId = clutch.clutch_id;
    const records = recordsByClutch[clutchId];
    for (const stage of STAGE_ORDER) {
      const stageCount = clutch.stage_counts[stage];
      if (
        stageCount.status !== "observed" ||
        stageCount.coverage !== "complete"
      ) {
        continue;
      }
      let derived;
      if (stage === "eggs_laid") {
        derived = records.length;
      } else {
        const [statusField, countedStatus] = statusFields[stage];
        for (const record of records) {
          if (record[statusField] === "unknown") {
            fail(
              `offspring_stage_status_unknown:${clutchId}:${stage}:` +
                record.offspring_record_id
            );
          }
        }
        derived = records.filter(
          (record) => record[statusField] === countedStatus
        ).length;
      }
      const declared = stageCount.count;
      if (declared !== derived) {
        fail(
          `stage_offspring_count_conflict:${clutchId}:${stage}:` +
            `${declared}:${derived}`
        );
      }
    }
  }
}

export function normalizeBreedingRecords(request) {
  const item = mapping(request, "breeding_record");
  rejectUnknownFields(item, TOP_LEVEL_FIELDS, "breeding_record");
  const schemaVersion = item.schema_version;
  if (schemaVersion !== "1.0.0") {
    fail(`breeding_schema_version_unsupported:${schemaVersion}`);
  }
  const caseId = requiredText(item.case_id, "case_id");
  const missingData = [];
  const { pairings, pairingIds } = normalizePairings(item.pairings);
  const { clutches, clutchIds } = normalizeClutches(
    item.clutches,
    pairingIds,
    missingData
  );
  const offspring = normalizeOffspring(
    item.offspring,
    clutchIds,
    missingData
  );
  reconcileCompleteCounts(clutches, offspring);
  return {
    schema_version: schemaVersion,
    case_id: caseId,
    pairings,
    clutches,
    offspring,
    missing_data: missingData,
    evidence: {
      source: "submitted_breeding_records",
      genetics_grade: "not_assessed",
      record_scope: "normalized_observations_only"
    },
    limits: {
      purpose: "breeding_record_normalization_only",
      gene_status: "not_assessed",
      inheritance_status: "not_assessed",
      decision_support_only: true,
      sex_method_recorded_not_verified: true
    }
  };
}
