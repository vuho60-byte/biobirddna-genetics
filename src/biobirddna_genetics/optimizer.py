from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Mapping, Sequence

from .cross import cross_uncertain
from .enums import Sex
from .errors import ValidationError
from .genotype import normalize_distribution
from .inference import _neumaier_sum, _portable_log2
from .models import GeneticProfile, ModelHypothesis, PairingWeights
from .pedigree import Pedigree


def _entropy(probabilities: Mapping[str, float]) -> float:
    # Portable log2 + sorted Neumaier sum: bit-identical with the Node
    # rankEntropy implementation on every platform and Python version.
    terms = sorted(
        probability * _portable_log2(probability)
        for probability in probabilities.values()
        if probability > 0.0
    )
    return -_neumaier_sum(terms)


def _normalize_hypotheses(
    hypotheses: Sequence[ModelHypothesis],
) -> list[ModelHypothesis]:
    if not hypotheses:
        raise ValidationError("model_hypotheses_required")
    total = sum(hypothesis.weight for hypothesis in hypotheses)
    if total <= 0.0:
        raise ValidationError("model_hypothesis_weights_zero")
    return [
        ModelHypothesis(model=hypothesis.model, weight=hypothesis.weight / total)
        for hypothesis in hypotheses
    ]


def _genotype_uncertainty(distribution: Mapping[str, float]) -> float:
    normalized = normalize_distribution(dict(distribution))
    if len(normalized) <= 1:
        return 0.0
    return _entropy(normalized) / _portable_log2(float(len(normalized)))


def _model_predictions_for_pair(
    male: GeneticProfile,
    female: GeneticProfile,
    hypotheses: Sequence[ModelHypothesis],
) -> tuple[list[tuple[float, dict[str, float]]], float, float]:
    predictions: list[tuple[float, dict[str, float]]] = []
    weighted_target_probability = 0.0
    weighted_genotype_uncertainty = 0.0

    for hypothesis in hypotheses:
        model = hypothesis.model
        male_distribution = male.genotype_posteriors.get(model.model_id)
        female_distribution = female.genotype_posteriors.get(model.model_id)
        if not male_distribution or not female_distribution:
            continue
        result = cross_uncertain(model, male_distribution, female_distribution)
        phenotype_probabilities = result["phenotype_probabilities"]
        predictions.append((hypothesis.weight, phenotype_probabilities))
        target_probability = sum(
            probability
            for key, probability in phenotype_probabilities.items()
            if key.endswith(f":{model.target_label}")
        )
        weighted_target_probability += hypothesis.weight * target_probability
        weighted_genotype_uncertainty += hypothesis.weight * 0.5 * (
            _genotype_uncertainty(male_distribution)
            + _genotype_uncertainty(female_distribution)
        )

    total_weight = sum(weight for weight, _ in predictions)
    if total_weight <= 0.0:
        return [], 0.0, 1.0
    predictions = [
        (weight / total_weight, distribution) for weight, distribution in predictions
    ]
    return (
        predictions,
        weighted_target_probability / total_weight,
        weighted_genotype_uncertainty / total_weight,
    )


def expected_information_gain(
    model_predictions: Sequence[tuple[float, Mapping[str, float]]],
) -> float:
    """Mutual information between the candidate model and one offspring outcome."""

    if not model_predictions:
        return 0.0
    categories = sorted(
        {
            category
            for _weight, probabilities in model_predictions
            for category in probabilities
        }
    )
    mixture: dict[str, float] = defaultdict(float)
    conditional_entropy = 0.0
    for weight, probabilities in model_predictions:
        normalized = {category: float(probabilities.get(category, 0.0)) for category in categories}
        total = sum(normalized.values())
        if total <= 0.0:
            continue
        normalized = {key: value / total for key, value in normalized.items()}
        for category, probability in normalized.items():
            mixture[category] += weight * probability
        conditional_entropy += weight * _entropy(normalized)
    return max(0.0, _entropy(mixture) - conditional_entropy)


def rank_mating_pairs(
    males: Sequence[GeneticProfile],
    females: Sequence[GeneticProfile],
    hypotheses: Sequence[ModelHypothesis],
    pedigree: Pedigree,
    *,
    weights: PairingWeights | None = None,
    max_offspring_inbreeding: float = 0.0625,
    block_direct_ancestor_and_full_sibling: bool = True,
) -> dict[str, Any]:
    weights = weights or PairingWeights()
    hypotheses = _normalize_hypotheses(hypotheses)
    if not 0.0 <= max_offspring_inbreeding <= 1.0:
        raise ValidationError("max_offspring_inbreeding_invalid")

    results: list[dict[str, Any]] = []
    for male in males:
        if male.sex != Sex.MALE:
            raise ValidationError(f"profile_not_male:{male.individual_id}")
        for female in females:
            if female.sex != Sex.FEMALE:
                raise ValidationError(f"profile_not_female:{female.individual_id}")

            block_reasons: list[str] = []
            if male.welfare_blocked or female.welfare_blocked:
                block_reasons.append("welfare_blocked")
            if male.individual_id == female.individual_id:
                block_reasons.append("same_individual_pair")

            hard_relation = pedigree.hard_relation_block(
                male.individual_id,
                female.individual_id,
            )
            if block_direct_ancestor_and_full_sibling and hard_relation:
                block_reasons.append(hard_relation)

            offspring_inbreeding = pedigree.offspring_inbreeding(
                male.individual_id,
                female.individual_id,
            )
            if offspring_inbreeding > max_offspring_inbreeding:
                block_reasons.append("offspring_inbreeding_above_configured_limit")

            predictions, target_probability, uncertainty = _model_predictions_for_pair(
                male,
                female,
                hypotheses,
            )
            if not predictions:
                block_reasons.append("missing_genotype_posteriors")
            information_gain = expected_information_gain(predictions)
            health_risk = max(
                male.health_risk,
                female.health_risk,
                male.fertility_risk,
                female.fertility_risk,
            )
            score = (
                weights.target * target_probability
                + weights.information * information_gain
                - weights.inbreeding * offspring_inbreeding
                - weights.health * health_risk
                - weights.uncertainty * uncertainty
            )
            if not math.isfinite(score):
                raise ValidationError(
                    f"pair_score_non_finite:{male.individual_id}:{female.individual_id}"
                )
            results.append(
                {
                    "male_id": male.individual_id,
                    "female_id": female.individual_id,
                    "eligible": not block_reasons,
                    "block_reasons": block_reasons,
                    "score": score,
                    "target_probability": target_probability,
                    "expected_information_gain_bits": information_gain,
                    "offspring_inbreeding": offspring_inbreeding,
                    "health_risk": health_risk,
                    "genotype_uncertainty": uncertainty,
                }
            )

    results.sort(
        key=lambda result: (
            result["eligible"],
            result["score"],
        ),
        reverse=True,
    )
    return {
        "pair_ranking": results,
        "configuration": {
            "max_offspring_inbreeding": max_offspring_inbreeding,
            "weights": {
                "target": weights.target,
                "information": weights.information,
                "inbreeding": weights.inbreeding,
                "health": weights.health,
                "uncertainty": weights.uncertainty,
            },
        },
        "limits": {
            "purpose": "decision_support_not_automatic_breeding_order",
            "health_and_welfare_review_required": True,
            "gene_status": "model_based",
        },
    }
