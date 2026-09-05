from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

from .breeding_records import normalize_breeding_records
from .cross import cross_known, cross_uncertain
from .enums import EvidenceGrade, Sex
from .errors import EvidenceBoundaryError, ValidationError
from .evidence import assert_claim_allowed, safe_claim_language
from .genotype import normalize_distribution
from .inference import default_candidate_models, infer_candidate_models
from .models import (
    GeneticProfile,
    InheritanceModel,
    ModelHypothesis,
    OffspringObservation,
    PairingWeights,
    PedigreeIndividual,
    PigmentObservation,
)
from .optimizer import rank_mating_pairs
from .pedigree import Pedigree
from .pigment import rank_pigment_hypotheses
from .simulation import simulate_recovery

INFERENCE_ROOT_FIELDS = {
    "target_label",
    "male_parent_phenotype",
    "female_parent_phenotype",
    "models",
    "model_priors",
    "offspring",
}
INFERENCE_OBSERVATION_FIELDS = {"sex", "phenotype", "count"}
PIGMENT_FIELDS = {
    "white_area_ratio",
    "brown_shift_score",
    "grey_shift_score",
    "dark_expansion_score",
    "red_patch_change_score",
    "pattern_retained",
    "iris_red_confirmed",
    "flash_red_eye_risk",
    "progressive_across_molts",
    "present_in_first_plumage",
    "localized_after_injury",
    "bilateral_symmetry",
    "image_quality",
}
CROSS_ROOT_FIELDS = {
    "model",
    "male_genotype",
    "female_genotype",
    "male_genotypes",
    "female_genotypes",
}
RANK_ROOT_FIELDS = {
    "pedigree",
    "males",
    "females",
    "hypotheses",
    "max_offspring_inbreeding",
    "block_direct_ancestor_and_full_sibling",
    "weights",
}
PEDIGREE_FIELDS = {"individual_id", "sex", "sire_id", "dam_id"}
PROFILE_FIELDS = {
    "individual_id",
    "sex",
    "genotype_posteriors",
    "phenotype_label",
    "health_risk",
    "fertility_risk",
    "welfare_blocked",
    "notes",
}
HYPOTHESIS_FIELDS = {"model", "weight"}
PAIRING_WEIGHT_FIELDS = {
    "target",
    "information",
    "inbreeding",
    "health",
    "uncertainty",
}
MAX_OFFSPRING_OBSERVATION_COUNT = 1_000_000

# Strongest claim each public operation may emit, paired with the evidence
# grade its output can maximally support (references/scientific-boundaries.md
# sections 1 and 4). Grades are the existing EvidenceGrade G0-G6; claim types
# are the vocabulary in evidence._MINIMUM. _apply_claim_boundary() runs
# assert_claim_allowed on the emitted language, so remapping an operation to a
# claim above its grade raises EvidenceBoundaryError before any result leaves
# dispatch() — the boundary is enforced by code, not by wrapper discipline.
# Mirrored in npm/src/index.mjs (OPERATION_CLAIM_CEILINGS) — keep in sync.
_OPERATION_CLAIM_CEILINGS: dict[str, tuple[str, EvidenceGrade]] = {
    # Pure simulation under assumed models: no genetic evidence at all.
    "cross": ("simulate_cross", EvidenceGrade.G0),
    # One submitted offspring dataset can at most support a candidate model
    # (G2 = repeated broods from the same pair; independence of pairs is not
    # verifiable from a single request).
    "infer": ("candidate_inheritance_model", EvidenceGrade.G2),
    # Visible phenotype only (G0 = image/appearance only).
    "pigment": ("phenotype_hypothesis", EvidenceGrade.G0),
    # Decision support from simulated crosses over hypothesized models.
    "rank-mates": ("simulate_cross", EvidenceGrade.G0),
    # Record normalization: observations recorded, nothing assessed.
    "normalize-records": ("observation_record", EvidenceGrade.G0),
    # Synthetic simulation only (output already carries evidence_grade G0).
    "simulate-recovery": ("simulate_cross", EvidenceGrade.G0),
}


def _apply_claim_boundary(
    operation: str,
    result: dict[str, Any],
) -> dict[str, Any]:
    ceiling = _OPERATION_CLAIM_CEILINGS.get(operation)
    if ceiling is None:
        raise EvidenceBoundaryError(
            f"operation_claim_ceiling_missing:{operation}"
        )
    claim_type, grade = ceiling
    language = safe_claim_language(claim_type, grade)
    assert_claim_allowed(language, grade)
    limits = result.get("limits")
    if not isinstance(limits, dict):
        raise EvidenceBoundaryError(f"result_limits_missing:{operation}")
    limits["claim_language"] = language
    return result


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValidationError(f"{name}_must_be_object")
    return value


def _reject_unknown_fields(
    value: Mapping[str, Any],
    allowed: set[str],
    name: str,
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValidationError(f"{name}_unknown_field:{unknown[0]}")


def _required_exact_text(value: Any, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
    ):
        raise ValidationError(f"{name}_required")
    return value


def _sequence(value: Any, name: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ValidationError(f"{name}_must_be_array")
    return value


def _reject_non_finite_numbers(
    value: Any,
    seen: set[int] | None = None,
) -> None:
    seen = seen if seen is not None else set()
    if isinstance(value, (Mapping, list, tuple)):
        identity = id(value)
        if identity in seen:
            raise ValidationError("input_json_cycle_forbidden")
        seen.add(identity)
        children = value.values() if isinstance(value, Mapping) else value
        for child in children:
            _reject_non_finite_numbers(child, seen)
        seen.remove(identity)
    elif (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and not math.isfinite(float(value))
    ):
        raise ValidationError("input_json_non_finite:number")


def handle_cross_request(request: Mapping[str, Any]) -> dict[str, Any]:
    _reject_unknown_fields(request, CROSS_ROOT_FIELDS, "cross_request")
    model = InheritanceModel.from_dict(_mapping(request.get("model"), "model"))
    male = request.get("male_genotype")
    female = request.get("female_genotype")
    male_distribution = request.get("male_genotypes")
    female_distribution = request.get("female_genotypes")

    known_mode = male is not None or female is not None
    uncertain_mode = (
        male_distribution is not None or female_distribution is not None
    )
    if known_mode and uncertain_mode:
        raise ValidationError("cross_parent_input_modes_mutually_exclusive")
    if male_distribution is not None or female_distribution is not None:
        if male_distribution is None or female_distribution is None:
            raise ValidationError("both_parent_genotype_distributions_required")
        return cross_uncertain(
            model,
            _mapping(male_distribution, "male_genotypes"),
            _mapping(female_distribution, "female_genotypes"),
        )
    if not isinstance(male, str) or not isinstance(female, str):
        raise ValidationError("male_and_female_genotype_required")
    return cross_known(model, male, female)


def handle_inference_request(request: Mapping[str, Any]) -> dict[str, Any]:
    _reject_unknown_fields(
        request,
        INFERENCE_ROOT_FIELDS,
        "inference_request",
    )
    target_label = _required_exact_text(
        request.get("target_label", "target_colour"),
        "target_label",
    )
    male_parent_phenotype = _required_exact_text(
        request.get("male_parent_phenotype"),
        "male_parent_phenotype",
    )
    female_parent_phenotype = _required_exact_text(
        request.get("female_parent_phenotype"),
        "female_parent_phenotype",
    )
    raw_models = request.get("models")
    models = (
        [
            InheritanceModel.from_dict(_mapping(model, "candidate_model"))
            for model in _sequence(raw_models, "models")
        ]
        if raw_models is not None
        else default_candidate_models(target_label)
    )
    observations: list[OffspringObservation] = []
    for raw_item in _sequence(request.get("offspring"), "offspring"):
        item = _mapping(raw_item, "offspring_observation")
        _reject_unknown_fields(
            item,
            INFERENCE_OBSERVATION_FIELDS,
            "offspring_observation",
        )
        phenotype = _required_exact_text(
            item.get("phenotype"),
            "offspring_phenotype",
        )
        count = item.get("count")
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or not 0 <= count <= MAX_OFFSPRING_OBSERVATION_COUNT
        ):
            raise ValidationError("offspring_count_invalid")
        raw_sex = item.get("sex", "unknown")
        if not isinstance(raw_sex, str):
            raise ValidationError("offspring_sex_invalid")
        try:
            sex = Sex(raw_sex)
        except ValueError as exc:
            raise ValidationError("offspring_sex_invalid") from exc
        observations.append(
            OffspringObservation(
                phenotype=phenotype,
                count=count,
                sex=sex,
            )
        )
    return infer_candidate_models(
        models,
        male_parent_phenotype=male_parent_phenotype,
        female_parent_phenotype=female_parent_phenotype,
        observations=observations,
        model_priors=(
            _mapping(request["model_priors"], "model_priors")
            if "model_priors" in request
            else None
        ),
    )


def handle_pigment_request(request: Mapping[str, Any]) -> dict[str, Any]:
    _reject_unknown_fields(request, PIGMENT_FIELDS, "pigment_request")

    def probability(name: str, default: float = 0.0) -> float:
        value = request.get(name, default)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or not 0.0 <= float(value) <= 1.0
        ):
            raise ValidationError(f"{name}_invalid")
        return float(value)

    def boolean(name: str, default: bool = False) -> bool:
        value = request.get(name, default)
        if not isinstance(value, bool):
            raise ValidationError(f"{name}_invalid")
        return value

    def nullable_boolean(name: str) -> bool | None:
        value = request.get(name)
        if value is not None and not isinstance(value, bool):
            raise ValidationError(f"{name}_invalid")
        return value

    image_quality = request.get("image_quality", "Q0")
    if image_quality not in {"Q0", "Q1", "Q2", "Q3", "Q4"}:
        raise ValidationError("image_quality_invalid")

    observation = PigmentObservation(
        white_area_ratio=probability("white_area_ratio"),
        brown_shift_score=probability("brown_shift_score"),
        grey_shift_score=probability("grey_shift_score"),
        dark_expansion_score=probability("dark_expansion_score"),
        red_patch_change_score=probability("red_patch_change_score"),
        pattern_retained=nullable_boolean("pattern_retained"),
        iris_red_confirmed=boolean("iris_red_confirmed"),
        flash_red_eye_risk=boolean("flash_red_eye_risk"),
        progressive_across_molts=boolean("progressive_across_molts"),
        present_in_first_plumage=nullable_boolean(
            "present_in_first_plumage"
        ),
        localized_after_injury=boolean("localized_after_injury"),
        bilateral_symmetry=(
            probability("bilateral_symmetry")
            if request.get("bilateral_symmetry") is not None
            else None
        ),
        image_quality=image_quality,
    )
    return rank_pigment_hypotheses(observation)


def _profile_from_dict(value: Mapping[str, Any]) -> GeneticProfile:
    _reject_unknown_fields(value, PROFILE_FIELDS, "genetic_profile")
    individual_id = _required_exact_text(
        value.get("individual_id"),
        "individual_id",
    )
    raw_sex = value.get("sex")
    if not isinstance(raw_sex, str):
        raise ValidationError("profile_sex_invalid")
    try:
        sex = Sex(raw_sex)
    except ValueError as exc:
        raise ValidationError("profile_sex_invalid") from exc
    raw_posteriors = _mapping(
        value.get("genotype_posteriors"),
        "genotype_posteriors",
    )

    def probability(name: str) -> float:
        raw = value.get(name, 0.0)
        if (
            isinstance(raw, bool)
            or not isinstance(raw, (int, float))
            or not math.isfinite(float(raw))
            or not 0.0 <= float(raw) <= 1.0
        ):
            raise ValidationError(f"{name}_invalid")
        return float(raw)

    welfare_blocked = value.get("welfare_blocked", False)
    if not isinstance(welfare_blocked, bool):
        raise ValidationError("welfare_blocked_invalid")
    phenotype_label = value.get("phenotype_label", "unknown")
    if not isinstance(phenotype_label, str):
        raise ValidationError("phenotype_label_invalid")
    raw_notes = value.get("notes", [])
    notes = _sequence(raw_notes, "notes")
    if any(not isinstance(item, str) for item in notes):
        raise ValidationError("notes_invalid")

    return GeneticProfile(
        individual_id=individual_id,
        sex=sex,
        genotype_posteriors={
            _required_exact_text(model_id, "model_id"): normalize_distribution(
                dict(
                    _mapping(
                    distribution,
                    "genotype_posterior",
                    )
                )
            )
            for model_id, distribution in raw_posteriors.items()
        },
        phenotype_label=phenotype_label,
        health_risk=probability("health_risk"),
        fertility_risk=probability("fertility_risk"),
        welfare_blocked=welfare_blocked,
        notes=list(notes),
    )


def handle_mate_ranking_request(request: Mapping[str, Any]) -> dict[str, Any]:
    _reject_unknown_fields(request, RANK_ROOT_FIELDS, "rank_mates_request")
    raw_maximum = request.get("max_offspring_inbreeding", 0.0625)
    if (
        isinstance(raw_maximum, bool)
        or not isinstance(raw_maximum, (int, float))
        or not math.isfinite(float(raw_maximum))
        or not 0.0 <= float(raw_maximum) <= 1.0
    ):
        raise ValidationError("max_offspring_inbreeding_invalid")
    raw_hard_block = request.get(
        "block_direct_ancestor_and_full_sibling",
        True,
    )
    if not isinstance(raw_hard_block, bool):
        raise ValidationError(
            "block_direct_ancestor_and_full_sibling_invalid"
        )
    raw_weights = (
        _mapping(request["weights"], "weights")
        if "weights" in request
        else None
    )
    if raw_weights is not None:
        _reject_unknown_fields(raw_weights, PAIRING_WEIGHT_FIELDS, "weights")
    pedigree_individuals: list[PedigreeIndividual] = []
    for raw_item in _sequence(request.get("pedigree"), "pedigree"):
        item = _mapping(raw_item, "pedigree_individual")
        _reject_unknown_fields(item, PEDIGREE_FIELDS, "pedigree_individual")
        individual_id = _required_exact_text(
            item.get("individual_id"),
            "individual_id",
        )
        raw_sex = item.get("sex", "unknown")
        if not isinstance(raw_sex, str):
            raise ValidationError("pedigree_sex_invalid")
        try:
            sex = Sex(raw_sex)
        except ValueError as exc:
            raise ValidationError("pedigree_sex_invalid") from exc

        def parent_id(name: str) -> str | None:
            raw_parent = item.get(name)
            if raw_parent is None:
                return None
            return _required_exact_text(raw_parent, name)

        pedigree_individuals.append(
            PedigreeIndividual(
                individual_id=individual_id,
                sex=sex,
                sire_id=parent_id("sire_id"),
                dam_id=parent_id("dam_id"),
            )
        )
    pedigree = Pedigree(pedigree_individuals)
    males = [
        _profile_from_dict(_mapping(item, "male_profile"))
        for item in _sequence(request.get("males"), "males")
    ]
    females = [
        _profile_from_dict(_mapping(item, "female_profile"))
        for item in _sequence(request.get("females"), "females")
    ]
    hypotheses: list[ModelHypothesis] = []
    for raw_item in _sequence(request.get("hypotheses"), "hypotheses"):
        item = _mapping(raw_item, "model_hypothesis")
        _reject_unknown_fields(
            item,
            HYPOTHESIS_FIELDS,
            "model_hypothesis",
        )
        if "weight" not in item:
            raise ValidationError("model_weight_invalid")
        hypotheses.append(
            ModelHypothesis(
                model=InheritanceModel.from_dict(
                    _mapping(item.get("model"), "model")
                ),
                weight=item["weight"],
            )
        )
    return rank_mating_pairs(
        males,
        females,
        hypotheses,
        pedigree,
        weights=PairingWeights.from_dict(raw_weights),
        max_offspring_inbreeding=float(raw_maximum),
        block_direct_ancestor_and_full_sibling=raw_hard_block,
    )


def handle_breeding_record_request(request: Mapping[str, Any]) -> dict[str, Any]:
    return normalize_breeding_records(request)


def handle_simulation_recovery_request(
    request: Mapping[str, Any],
) -> dict[str, Any]:
    return simulate_recovery(request)


def dispatch(operation: str, request: Mapping[str, Any]) -> dict[str, Any]:
    _reject_non_finite_numbers(request)
    if operation == "cross":
        result = handle_cross_request(request)
    elif operation == "infer":
        result = handle_inference_request(request)
    elif operation == "pigment":
        result = handle_pigment_request(request)
    elif operation == "rank-mates":
        result = handle_mate_ranking_request(request)
    elif operation == "normalize-records":
        result = handle_breeding_record_request(request)
    elif operation == "simulate-recovery":
        result = handle_simulation_recovery_request(request)
    else:
        raise ValidationError(f"operation_unknown:{operation}")
    return _apply_claim_boundary(operation, result)
