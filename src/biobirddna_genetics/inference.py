from __future__ import annotations

import math
from collections import defaultdict
from itertools import product
from typing import Any, Iterable, Mapping, Sequence

from .cross import collapse_phenotypes_by_sex, cross_known
from .enums import InheritanceMode, Sex
from .errors import ValidationError
from .genotype import compatible_genotypes
from .models import InheritanceModel, OffspringObservation


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        return float("-inf")
    maximum = max(values)
    if math.isinf(maximum):
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def _normalize_log_weights(values: Mapping[str, float]) -> dict[str, float]:
    denominator = _logsumexp(list(values.values()))
    if math.isinf(denominator):
        uniform = 1.0 / len(values)
        return {key: uniform for key in values}
    return {key: math.exp(value - denominator) for key, value in values.items()}


def _json_safe_numbers(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _json_safe_numbers(child)
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [_json_safe_numbers(child) for child in value]
    if (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and not math.isfinite(value)
    ):
        return None
    return value


def _observation_mode(observations: Iterable[OffspringObservation]) -> str:
    sexes = {observation.sex for observation in observations if observation.count > 0}
    if not sexes or sexes == {Sex.UNKNOWN}:
        return "sex_collapsed"
    if Sex.UNKNOWN in sexes:
        raise ValidationError("cannot_mix_known_and_unknown_offspring_sex")
    return "sex_specific"


def _counts(observations: Iterable[OffspringObservation], mode: str) -> dict[str, int]:
    result: dict[str, int] = defaultdict(int)
    for observation in observations:
        key = (
            observation.phenotype
            if mode == "sex_collapsed"
            else f"{observation.sex.value}:{observation.phenotype}"
        )
        result[key] += observation.count
    return dict(result)


def _smooth_probabilities(
    probabilities: Mapping[str, float],
    categories: Sequence[str],
    error_rate: float,
) -> dict[str, float]:
    if not categories:
        raise ValidationError("observation_categories_empty")
    k = len(categories)
    smoothed = {
        category: (1.0 - error_rate) * float(probabilities.get(category, 0.0))
        + error_rate / k
        for category in categories
    }
    total = sum(smoothed.values())
    if total <= 0.0:
        raise ValidationError("smoothed_probability_zero_total")
    return {key: value / total for key, value in smoothed.items()}


def multinomial_log_likelihood(
    counts: Mapping[str, int],
    probabilities: Mapping[str, float],
) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    coefficient = math.lgamma(total + 1.0) - sum(
        math.lgamma(count + 1.0) for count in counts.values()
    )
    likelihood = coefficient
    for category, count in counts.items():
        probability = float(probabilities.get(category, 0.0))
        if count > 0 and probability <= 0.0:
            return float("-inf")
        if count > 0:
            likelihood += count * math.log(probability)
    return likelihood


def _configuration_predictions(
    model: InheritanceModel,
    male_genotype: str,
    female_genotype: str,
    mode: str,
) -> dict[str, float]:
    result = cross_known(model, male_genotype, female_genotype)
    probabilities = result["phenotype_probabilities"]
    return (
        collapse_phenotypes_by_sex(probabilities)
        if mode == "sex_collapsed"
        else probabilities
    )


def _score_non_genetic(
    model: InheritanceModel,
    counts: Mapping[str, int],
    mode: str,
) -> tuple[float, list[dict[str, Any]]]:
    p = model.baseline_target_probability
    if mode == "sex_collapsed":
        raw = {model.target_label: p, model.normal_label: 1.0 - p}
    else:
        female = model.female_probability
        raw = {
            f"male:{model.target_label}": (1.0 - female) * p,
            f"male:{model.normal_label}": (1.0 - female) * (1.0 - p),
            f"female:{model.target_label}": female * p,
            f"female:{model.normal_label}": female * (1.0 - p),
        }
    categories = sorted(set(counts) | set(raw))
    probabilities = _smooth_probabilities(raw, categories, model.observation_error)
    log_likelihood = multinomial_log_likelihood(counts, probabilities)
    return log_likelihood, [
        {
            "male_genotype": "NA",
            "female_genotype": "NA",
            "prior": 1.0,
            "posterior_within_model": 1.0,
            "log_likelihood": log_likelihood,
            "predicted_probabilities": probabilities,
        }
    ]


def score_model(
    model: InheritanceModel,
    male_parent_phenotype: str,
    female_parent_phenotype: str,
    observations: Sequence[OffspringObservation],
    *,
    male_genotype_prior: Mapping[str, float] | None = None,
    female_genotype_prior: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    mode = _observation_mode(observations)
    observed_counts = _counts(observations, mode)

    if model.mode == InheritanceMode.NON_GENETIC:
        log_evidence, configurations = _score_non_genetic(
            model,
            observed_counts,
            mode,
        )
        return {
            "model_id": model.model_id,
            "mode": model.mode.value,
            "log_evidence": log_evidence,
            "observation_mode": mode,
            "observed_counts": observed_counts,
            "configurations": configurations,
            "best_parent_configuration": configurations[0],
        }

    male_genotypes = compatible_genotypes(
        model,
        Sex.MALE,
        male_parent_phenotype,
    )
    female_genotypes = compatible_genotypes(
        model,
        Sex.FEMALE,
        female_parent_phenotype,
    )
    if not male_genotypes or not female_genotypes:
        return {
            "model_id": model.model_id,
            "mode": model.mode.value,
            "log_evidence": float("-inf"),
            "observation_mode": mode,
            "observed_counts": observed_counts,
            "configurations": [],
            "best_parent_configuration": None,
            "reason": "parent_phenotype_incompatible_with_model",
        }

    def prior_for(
        genotype: str,
        allowed: Sequence[str],
        explicit: Mapping[str, float] | None,
    ) -> float:
        if explicit is None:
            return 1.0 / len(allowed)
        raw = float(explicit.get(genotype, 0.0))
        total = sum(float(explicit.get(item, 0.0)) for item in allowed)
        return raw / total if total > 0.0 else 0.0

    configurations: list[dict[str, Any]] = []
    log_terms: list[float] = []
    for male_genotype, female_genotype in product(
        male_genotypes,
        female_genotypes,
    ):
        male_prior = prior_for(male_genotype, male_genotypes, male_genotype_prior)
        female_prior = prior_for(
            female_genotype,
            female_genotypes,
            female_genotype_prior,
        )
        prior = male_prior * female_prior
        if prior <= 0.0:
            continue
        raw_probabilities = _configuration_predictions(
            model,
            male_genotype,
            female_genotype,
            mode,
        )
        categories = sorted(set(observed_counts) | set(raw_probabilities))
        probabilities = _smooth_probabilities(
            raw_probabilities,
            categories,
            model.observation_error,
        )
        log_likelihood = multinomial_log_likelihood(
            observed_counts,
            probabilities,
        )
        log_joint = math.log(prior) + log_likelihood
        log_terms.append(log_joint)
        configurations.append(
            {
                "male_genotype": male_genotype,
                "female_genotype": female_genotype,
                "prior": prior,
                "log_likelihood": log_likelihood,
                "log_joint": log_joint,
                "predicted_probabilities": probabilities,
            }
        )

    log_evidence = _logsumexp(log_terms)
    for configuration in configurations:
        configuration["posterior_within_model"] = (
            math.exp(configuration["log_joint"] - log_evidence)
            if not math.isinf(log_evidence)
            else 0.0
        )
    configurations.sort(
        key=lambda configuration: configuration["posterior_within_model"],
        reverse=True,
    )
    return {
        "model_id": model.model_id,
        "mode": model.mode.value,
        "log_evidence": log_evidence,
        "observation_mode": mode,
        "observed_counts": observed_counts,
        "configurations": configurations,
        "best_parent_configuration": configurations[0] if configurations else None,
    }


def infer_candidate_models(
    models: Sequence[InheritanceModel],
    male_parent_phenotype: str,
    female_parent_phenotype: str,
    observations: Sequence[OffspringObservation],
    *,
    model_priors: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    if not models:
        raise ValidationError("candidate_models_required")
    if not observations:
        raise ValidationError("offspring_observations_required")
    seen_model_ids: set[str] = set()
    for model in models:
        if model.model_id in seen_model_ids:
            raise ValidationError(
                f"candidate_model_id_duplicate:{model.model_id}"
            )
        seen_model_ids.add(model.model_id)
    m0_count = sum(
        model.model_id == "M0_NONGENETIC"
        and model.mode == InheritanceMode.NON_GENETIC
        for model in models
    )
    if m0_count != 1:
        raise ValidationError(
            "candidate_models_require_exactly_one_M0_NONGENETIC"
        )
    if len(models) < 2:
        raise ValidationError("candidate_models_require_multiple_models")

    raw_priors: dict[str, float] = {}
    for model in models:
        prior = (
            float(model_priors.get(model.model_id, 0.0))
            if model_priors is not None
            else 1.0
        )
        if not math.isfinite(prior):
            raise ValidationError(
                f"model_prior_non_finite:{model.model_id}"
            )
        if prior < 0.0:
            raise ValidationError(f"model_prior_negative:{model.model_id}")
        raw_priors[model.model_id] = prior
    prior_total = sum(raw_priors.values())
    if prior_total <= 0.0:
        raise ValidationError("model_priors_zero_total")
    priors = {key: value / prior_total for key, value in raw_priors.items()}

    scores: list[dict[str, Any]] = []
    log_posterior_unnormalized: dict[str, float] = {}
    for model in models:
        score = score_model(
            model,
            male_parent_phenotype,
            female_parent_phenotype,
            observations,
        )
        score["prior"] = priors[model.model_id]
        scores.append(score)
        log_posterior_unnormalized[model.model_id] = (
            score["log_evidence"] + math.log(priors[model.model_id])
            if priors[model.model_id] > 0.0
            else float("-inf")
        )

    posterior = _normalize_log_weights(log_posterior_unnormalized)
    for score in scores:
        score["posterior"] = posterior[score["model_id"]]
    scores.sort(
        key=lambda value: (-value["posterior"], value["model_id"])
    )

    return _json_safe_numbers({
        "ranking": scores,
        "best_model": scores[0]["model_id"],
        "best_model_posterior": scores[0]["posterior"],
        "limits": {
            "claim": "candidate_inheritance_model_only",
            "gene_status": "not_assessed",
            "minimum_language": "phu_hop_hon_voi_mo_hinh",
        },
    })


def default_candidate_models(
    target_label: str = "target_colour",
) -> list[InheritanceModel]:
    return [
        InheritanceModel(
            model_id="M1_AR",
            mode=InheritanceMode.AUTOSOMAL_RECESSIVE,
            target_label=target_label,
            description="One-locus autosomal recessive candidate",
        ),
        InheritanceModel(
            model_id="M2_AD",
            mode=InheritanceMode.AUTOSOMAL_DOMINANT,
            target_label=target_label,
            description="One-locus autosomal dominant candidate",
        ),
        InheritanceModel(
            model_id="M3_INC",
            mode=InheritanceMode.INCOMPLETE_DOMINANT,
            target_label=target_label,
            description="One-locus incomplete-dominance candidate",
        ),
        InheritanceModel(
            model_id="M4_ZR",
            mode=InheritanceMode.Z_LINKED_RECESSIVE,
            target_label=target_label,
            description="Z-linked recessive candidate under avian ZZ/ZW convention",
        ),
        InheritanceModel(
            model_id="M5_ZD",
            mode=InheritanceMode.Z_LINKED_DOMINANT,
            target_label=target_label,
            description="Z-linked dominant candidate under avian ZZ/ZW convention",
        ),
        InheritanceModel(
            model_id="M0_NONGENETIC",
            mode=InheritanceMode.NON_GENETIC,
            target_label=target_label,
            baseline_target_probability=0.05,
            description="Independent non-heritable baseline candidate",
        ),
    ]
