from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

from .errors import ValidationError

STAGE_ORDER = (
    "eggs_laid",
    "fertile_eggs",
    "hatched",
    "survived_to_first_molt",
)
STAGE_MISSING_IMPACT = {
    "eggs_laid": "egg_stage_not_observed",
    "fertile_eggs": "fertility_stage_not_observed",
    "hatched": "hatch_stage_not_observed",
    "survived_to_first_molt": "survival_stage_not_observed",
}
STAGE_SCOPE = {
    "eggs_laid": "egg",
    "fertile_eggs": "fertility",
    "hatched": "hatch",
    "survived_to_first_molt": "survival",
}
OUTCOMES = {
    "alive",
    "unfertilized",
    "embryo_death",
    "unhatched",
    "hatch_failure",
    "post_hatch_death",
    "lost_to_followup",
    "unknown",
}
SEXES = {"male", "female", "unknown"}
SEX_METHODS = {
    "dna",
    "reproductive",
    "veterinary",
    "visual",
    "reported",
    "unknown",
}
FERTILITY_STATUSES = {"fertile", "unfertilized", "unknown"}
HATCH_STATUSES = {"hatched", "not_hatched", "not_applicable", "unknown"}
SURVIVAL_STATUSES = {
    "survived",
    "died_before_first_molt",
    "not_applicable",
    "unknown",
}
VALID_STAGE_PATHS = {
    ("unfertilized", "not_applicable", "not_applicable"),
    ("fertile", "not_hatched", "not_applicable"),
    ("fertile", "hatched", "survived"),
    ("fertile", "hatched", "died_before_first_molt"),
    ("fertile", "hatched", "unknown"),
    ("fertile", "unknown", "unknown"),
    ("unknown", "unknown", "unknown"),
}
TOP_LEVEL_FIELDS = {
    "schema_version",
    "case_id",
    "pairings",
    "clutches",
    "offspring",
}
PAIRING_FIELDS = {
    "pairing_id",
    "male_parent_id",
    "female_parent_id",
    "male_parent_phenotype",
    "female_parent_phenotype",
    "male_parent_sex_method",
    "male_parent_sex_confidence",
    "female_parent_sex_method",
    "female_parent_sex_confidence",
    "reciprocal_group_id",
}
CLUTCH_FIELDS = {
    "clutch_id",
    "pairing_id",
    "offspring_records_coverage",
    "stage_counts",
}
STAGE_COUNT_FIELDS = {"status", "coverage", "count"}
PHENOTYPE_OBSERVATION_FIELDS = {
    "observation_id",
    "sequence_index",
    "stage",
    "status",
    "phenotype",
}
OFFSPRING_FIELDS = {
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
    "phenotype_timeline",
}


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValidationError(f"{name}_must_be_object")
    return value


def _sequence(value: Any, name: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ValidationError(f"{name}_must_be_array")
    return value


def _required_text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{name}_required")
    return value


def _reject_unknown_fields(
    item: Mapping[str, Any],
    allowed: set[str],
    name: str,
) -> None:
    unknown = sorted(set(item) - allowed)
    if unknown:
        raise ValidationError(f"{name}_unknown_field:{unknown[0]}")


def _optional_text(value: Any, name: str) -> str | None:
    if value is None:
        return None
    return _required_text(value, name)


def _unique_id(value: Any, name: str, seen: set[str]) -> str:
    identifier = _required_text(value, name)
    if identifier in seen:
        raise ValidationError(f"{name}_duplicate:{identifier}")
    seen.add(identifier)
    return identifier


def _normalize_stage_count(
    value: Any,
    *,
    clutch_id: str,
    stage: str,
    missing_data: list[dict[str, str]],
) -> dict[str, int | str | None]:
    item = _mapping(value, f"{stage}_stage_count")
    _reject_unknown_fields(item, STAGE_COUNT_FIELDS, "stage_count")
    status = item.get("status")
    if status == "missing":
        coverage = item.get("coverage")
        if coverage != "unknown":
            raise ValidationError(
                f"missing_stage_coverage_must_be_unknown:{clutch_id}:{stage}"
            )
        if item.get("count") is not None:
            raise ValidationError(
                f"missing_stage_count_must_be_null:{clutch_id}:{stage}"
            )
        missing_data.append(
            {
                "field": f"clutches/{clutch_id}/stage_counts/{stage}",
                "impact": STAGE_MISSING_IMPACT[stage],
            }
        )
        return {"status": "missing", "coverage": "unknown", "count": None}
    if status != "observed":
        raise ValidationError(f"stage_count_status_invalid:{clutch_id}:{stage}")
    coverage = item.get("coverage")
    if coverage not in {"complete", "partial"}:
        raise ValidationError(f"stage_coverage_invalid:{clutch_id}:{stage}")
    count = item.get("count")
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValidationError(f"stage_count_invalid:{clutch_id}:{stage}")
    if coverage == "partial":
        missing_data.append(
            {
                "field": f"clutches/{clutch_id}/stage_counts/{stage}",
                "impact": f"{STAGE_SCOPE[stage]}_stage_incomplete_coverage",
            }
        )
    return {"status": "observed", "coverage": coverage, "count": count}


def _normalize_pairings(value: Any) -> tuple[list[dict[str, Any]], set[str]]:
    pairings: list[dict[str, Any]] = []
    pairing_ids: set[str] = set()
    raw_pairings = _sequence(value, "pairings")
    if not raw_pairings:
        raise ValidationError("pairings_must_not_be_empty")
    for raw in raw_pairings:
        item = _mapping(raw, "pairing")
        _reject_unknown_fields(item, PAIRING_FIELDS, "pairing")
        pairing_id = _unique_id(item.get("pairing_id"), "pairing_id", pairing_ids)
        male_parent_id = _required_text(
            item.get("male_parent_id"),
            "male_parent_id",
        )
        female_parent_id = _required_text(
            item.get("female_parent_id"),
            "female_parent_id",
        )
        if male_parent_id == female_parent_id:
            raise ValidationError(f"pairing_parents_must_differ:{pairing_id}")
        parent_sex_provenance: dict[str, str | float] = {}
        for role in ("male", "female"):
            method = item.get(f"{role}_parent_sex_method")
            confidence_value = item.get(f"{role}_parent_sex_confidence")
            if (
                method not in SEX_METHODS - {"unknown"}
                or isinstance(confidence_value, bool)
                or not isinstance(confidence_value, (int, float))
                or not math.isfinite(float(confidence_value))
                or not 0.0 <= float(confidence_value) <= 1.0
            ):
                raise ValidationError(
                    f"parent_sex_provenance_required:{pairing_id}:{role}"
                )
            parent_sex_provenance[f"{role}_parent_sex_method"] = str(method)
            parent_sex_provenance[f"{role}_parent_sex_confidence"] = float(
                confidence_value
            )
        pairings.append(
            {
                "pairing_id": pairing_id,
                "male_parent_id": male_parent_id,
                "female_parent_id": female_parent_id,
                "male_parent_phenotype": _required_text(
                    item.get("male_parent_phenotype"),
                    "male_parent_phenotype",
                ),
                "female_parent_phenotype": _required_text(
                    item.get("female_parent_phenotype"),
                    "female_parent_phenotype",
                ),
                **parent_sex_provenance,
                "reciprocal_group_id": _optional_text(
                    item.get("reciprocal_group_id"),
                    "reciprocal_group_id",
                ),
                "cross_direction_key": f"{male_parent_id}>{female_parent_id}",
            }
        )
    pairings.sort(key=lambda item: item["pairing_id"])
    return pairings, pairing_ids


def _normalize_clutches(
    value: Any,
    pairing_ids: set[str],
    missing_data: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], set[str]]:
    clutches: list[dict[str, Any]] = []
    clutch_ids: set[str] = set()
    for raw in _sequence(value, "clutches"):
        item = _mapping(raw, "clutch")
        _reject_unknown_fields(item, CLUTCH_FIELDS, "clutch")
        clutch_id = _unique_id(item.get("clutch_id"), "clutch_id", clutch_ids)
        pairing_id = _required_text(item.get("pairing_id"), "pairing_id")
        if pairing_id not in pairing_ids:
            raise ValidationError(f"clutch_pairing_unknown:{clutch_id}:{pairing_id}")
        offspring_records_coverage = item.get("offspring_records_coverage")
        if offspring_records_coverage not in {"complete", "partial", "unknown"}:
            raise ValidationError(
                f"offspring_records_coverage_invalid:{clutch_id}"
            )
        if offspring_records_coverage != "complete":
            missing_data.append(
                {
                    "field": f"clutches/{clutch_id}/offspring",
                    "impact": (
                        "offspring_records_incomplete_coverage"
                        if offspring_records_coverage == "partial"
                        else "offspring_records_unknown_coverage"
                    ),
                }
            )
        raw_counts = _mapping(item.get("stage_counts"), "stage_counts")
        _reject_unknown_fields(raw_counts, set(STAGE_ORDER), "stage_counts")
        stage_counts: dict[str, dict[str, int | str | None]] = {}
        for stage in STAGE_ORDER:
            if stage not in raw_counts:
                raise ValidationError(f"stage_count_required:{clutch_id}:{stage}")
            stage_counts[stage] = _normalize_stage_count(
                raw_counts.get(stage),
                clutch_id=clutch_id,
                stage=stage,
                missing_data=missing_data,
            )
        previous_stage: str | None = None
        previous_count: int | None = None
        for stage in STAGE_ORDER:
            observation = stage_counts[stage]
            if (
                observation["status"] != "observed"
                or observation["coverage"] != "complete"
            ):
                continue
            count = int(observation["count"])
            if previous_count is not None and count > previous_count:
                raise ValidationError(
                    "clutch_stage_count_increases:"
                    f"{clutch_id}:{previous_stage}:{previous_count}:{stage}:{count}"
                )
            previous_stage = stage
            previous_count = count
        clutches.append(
            {
                "clutch_id": clutch_id,
                "pairing_id": pairing_id,
                "offspring_records_coverage": offspring_records_coverage,
                "stage_counts": stage_counts,
            }
        )
    clutches.sort(key=lambda item: item["clutch_id"])
    return clutches, clutch_ids


def _confidence(value: Any, offspring_record_id: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValidationError(f"sex_confidence_invalid:{offspring_record_id}")
    confidence = float(value)
    if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
        raise ValidationError(f"sex_confidence_invalid:{offspring_record_id}")
    return confidence


def _normalize_timeline(
    value: Any,
    offspring_record_id: str,
    seen_observation_ids: set[str],
    missing_data: list[dict[str, str]],
) -> list[dict[str, Any]]:
    timeline: list[dict[str, Any]] = []
    sequence_indices: set[int] = set()
    for raw in _sequence(value, "phenotype_timeline"):
        item = _mapping(raw, "phenotype_observation")
        _reject_unknown_fields(
            item,
            PHENOTYPE_OBSERVATION_FIELDS,
            "phenotype_observation",
        )
        observation_id = _unique_id(
            item.get("observation_id"),
            "phenotype_observation_id",
            seen_observation_ids,
        )
        sequence_index = item.get("sequence_index")
        if (
            isinstance(sequence_index, bool)
            or not isinstance(sequence_index, int)
            or sequence_index < 0
        ):
            raise ValidationError(
                f"phenotype_sequence_index_invalid:{observation_id}"
            )
        if sequence_index in sequence_indices:
            raise ValidationError(
                "phenotype_sequence_index_duplicate:"
                f"{offspring_record_id}:{sequence_index}"
            )
        sequence_indices.add(sequence_index)
        stage = _required_text(item.get("stage"), "phenotype_stage")
        status = item.get("status")
        if status == "observed":
            phenotype = _required_text(item.get("phenotype"), "phenotype")
        elif status == "missing":
            if item.get("phenotype") is not None:
                raise ValidationError(
                    f"missing_phenotype_must_be_null:{observation_id}"
                )
            phenotype = None
            missing_data.append(
                {
                    "field": (
                        f"offspring/{offspring_record_id}/"
                        f"phenotype_timeline/{observation_id}"
                    ),
                    "impact": "phenotype_stage_not_observed",
                }
            )
        else:
            raise ValidationError(f"phenotype_status_invalid:{observation_id}")
        timeline.append(
            {
                "observation_id": observation_id,
                "sequence_index": sequence_index,
                "stage": stage,
                "status": status,
                "phenotype": phenotype,
            }
        )
    timeline.sort(key=lambda item: item["sequence_index"])
    return timeline


def _stage_status(
    value: Any,
    *,
    allowed: set[str],
    name: str,
    offspring_record_id: str,
) -> str:
    if value not in allowed:
        raise ValidationError(f"{name}_invalid:{offspring_record_id}")
    return str(value)


def _validate_stage_path(
    *,
    offspring_record_id: str,
    outcome: str,
    fertility_status: str,
    hatch_status: str,
    survival_status: str,
) -> None:
    stage_path = (fertility_status, hatch_status, survival_status)
    if stage_path not in VALID_STAGE_PATHS:
        raise ValidationError(
            f"offspring_stage_path_invalid:{offspring_record_id}"
        )
    outcome_paths = {
        "alive": {
            ("fertile", "hatched", "survived"),
            ("fertile", "hatched", "unknown"),
        },
        "unfertilized": {
            ("unfertilized", "not_applicable", "not_applicable")
        },
        "embryo_death": {
            ("fertile", "not_hatched", "not_applicable")
        },
        "unhatched": {("fertile", "not_hatched", "not_applicable")},
        "hatch_failure": {("fertile", "not_hatched", "not_applicable")},
        "post_hatch_death": {
            ("fertile", "hatched", "died_before_first_molt")
        },
        "lost_to_followup": {
            ("fertile", "hatched", "unknown"),
            ("fertile", "unknown", "unknown"),
            ("unknown", "unknown", "unknown"),
        },
        "unknown": VALID_STAGE_PATHS,
    }
    if stage_path not in outcome_paths[outcome]:
        raise ValidationError(
            f"offspring_outcome_stage_conflict:{offspring_record_id}:{outcome}"
        )


def _normalize_offspring(
    value: Any,
    clutch_ids: set[str],
    missing_data: list[dict[str, str]],
) -> list[dict[str, Any]]:
    offspring: list[dict[str, Any]] = []
    offspring_ids: set[str] = set()
    observation_ids: set[str] = set()
    for raw in _sequence(value, "offspring"):
        item = _mapping(raw, "offspring_record")
        _reject_unknown_fields(item, OFFSPRING_FIELDS, "offspring_record")
        offspring_record_id = _unique_id(
            item.get("offspring_record_id"),
            "offspring_record_id",
            offspring_ids,
        )
        clutch_id = _required_text(item.get("clutch_id"), "clutch_id")
        if clutch_id not in clutch_ids:
            raise ValidationError(
                f"offspring_clutch_unknown:{offspring_record_id}:{clutch_id}"
            )
        outcome = item.get("outcome")
        if outcome not in OUTCOMES:
            raise ValidationError(f"offspring_outcome_invalid:{offspring_record_id}")
        for required_field in ("sex_confidence", "phenotype_timeline"):
            if required_field not in item:
                raise ValidationError(
                    "offspring_field_required:"
                    f"{offspring_record_id}:{required_field}"
                )
        fertility_status = _stage_status(
            item.get("fertility_status"),
            allowed=FERTILITY_STATUSES,
            name="fertility_status",
            offspring_record_id=offspring_record_id,
        )
        hatch_status = _stage_status(
            item.get("hatch_status"),
            allowed=HATCH_STATUSES,
            name="hatch_status",
            offspring_record_id=offspring_record_id,
        )
        survival_status = _stage_status(
            item.get("survival_status"),
            allowed=SURVIVAL_STATUSES,
            name="survival_status",
            offspring_record_id=offspring_record_id,
        )
        _validate_stage_path(
            offspring_record_id=offspring_record_id,
            outcome=str(outcome),
            fertility_status=fertility_status,
            hatch_status=hatch_status,
            survival_status=survival_status,
        )
        sex = item.get("sex")
        if sex not in SEXES:
            raise ValidationError(f"offspring_sex_invalid:{offspring_record_id}")
        sex_method = item.get("sex_method")
        if sex_method not in SEX_METHODS:
            raise ValidationError(f"sex_method_invalid:{offspring_record_id}")
        sex_confidence = _confidence(
            item.get("sex_confidence"),
            offspring_record_id,
        )
        if sex != "unknown" and (
            sex_method == "unknown" or sex_confidence is None
        ):
            raise ValidationError(
                f"known_sex_requires_method_and_confidence:{offspring_record_id}"
            )
        offspring.append(
            {
                "offspring_record_id": offspring_record_id,
                "individual_id": _optional_text(
                    item.get("individual_id"),
                    "individual_id",
                ),
                "clutch_id": clutch_id,
                "outcome": outcome,
                "fertility_status": fertility_status,
                "hatch_status": hatch_status,
                "survival_status": survival_status,
                "sex": sex,
                "sex_method": sex_method,
                "sex_confidence": sex_confidence,
                "phenotype_timeline": _normalize_timeline(
                    item.get("phenotype_timeline"),
                    offspring_record_id,
                    observation_ids,
                    missing_data,
                ),
            }
        )
    return offspring


def _reconcile_complete_counts(
    clutches: Sequence[Mapping[str, Any]],
    offspring: Sequence[Mapping[str, Any]],
) -> None:
    records_by_clutch: dict[str, list[Mapping[str, Any]]] = {
        str(clutch["clutch_id"]): [] for clutch in clutches
    }
    for record in offspring:
        records_by_clutch[str(record["clutch_id"])].append(record)

    status_fields = {
        "fertile_eggs": ("fertility_status", "fertile"),
        "hatched": ("hatch_status", "hatched"),
        "survived_to_first_molt": ("survival_status", "survived"),
    }
    for clutch in clutches:
        if clutch["offspring_records_coverage"] != "complete":
            continue
        clutch_id = str(clutch["clutch_id"])
        records = records_by_clutch[clutch_id]
        for stage in STAGE_ORDER:
            stage_count = clutch["stage_counts"][stage]
            if (
                stage_count["status"] != "observed"
                or stage_count["coverage"] != "complete"
            ):
                continue
            if stage == "eggs_laid":
                derived = len(records)
            else:
                status_field, counted_status = status_fields[stage]
                for record in records:
                    if record[status_field] == "unknown":
                        raise ValidationError(
                            "offspring_stage_status_unknown:"
                            f"{clutch_id}:{stage}:"
                            f"{record['offspring_record_id']}"
                        )
                derived = sum(
                    record[status_field] == counted_status for record in records
                )
            declared = int(stage_count["count"])
            if declared != derived:
                raise ValidationError(
                    "stage_offspring_count_conflict:"
                    f"{clutch_id}:{stage}:{declared}:{derived}"
                )


def normalize_breeding_records(request: Mapping[str, Any]) -> dict[str, Any]:
    request = _mapping(request, "breeding_record")
    _reject_unknown_fields(request, TOP_LEVEL_FIELDS, "breeding_record")
    schema_version = request.get("schema_version")
    if schema_version != "1.0.0":
        raise ValidationError(
            f"breeding_schema_version_unsupported:{schema_version}"
        )
    case_id = _required_text(request.get("case_id"), "case_id")
    missing_data: list[dict[str, str]] = []
    pairings, pairing_ids = _normalize_pairings(request.get("pairings"))
    clutches, clutch_ids = _normalize_clutches(
        request.get("clutches"),
        pairing_ids,
        missing_data,
    )
    offspring = _normalize_offspring(
        request.get("offspring"),
        clutch_ids,
        missing_data,
    )
    _reconcile_complete_counts(clutches, offspring)
    return {
        "schema_version": schema_version,
        "case_id": case_id,
        "pairings": pairings,
        "clutches": clutches,
        "offspring": offspring,
        "missing_data": missing_data,
        "evidence": {
            "source": "submitted_breeding_records",
            "genetics_grade": "not_assessed",
            "record_scope": "normalized_observations_only",
        },
        "limits": {
            "purpose": "breeding_record_normalization_only",
            "gene_status": "not_assessed",
            "inheritance_status": "not_assessed",
            "decision_support_only": True,
            "sex_method_recorded_not_verified": True,
        },
    }
