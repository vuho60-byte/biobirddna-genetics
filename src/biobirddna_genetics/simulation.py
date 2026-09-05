from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from typing import Any, Mapping, Sequence

from .cross import cross_known
from .enums import Sex
from .errors import ValidationError
from .inference import (
    _LN2,
    _neumaier_sum,
    _portable_log,
    _portable_log2,
    default_candidate_models,
    infer_candidate_models,
    score_model,
)
from .genotype import phenotype_distribution
from .models import InheritanceModel, OffspringObservation

ROOT_FIELDS = {
    "schema_version",
    "base_seed",
    "goodness_of_fit_replicates",
    "misspecification_alpha",
    "experiments",
}
EXPERIMENT_FIELDS = {
    "experiment_id",
    "generating_model",
    "male_genotype",
    "female_genotype",
    "male_parent_phenotype",
    "female_parent_phenotype",
    "brood_sizes",
    "replicates",
    "candidate_model_ids",
}
MODEL_FIELDS = {
    "model_id",
    "mode",
    "normal_allele",
    "variant_allele",
    "normal_label",
    "target_label",
    "intermediate_label",
    "penetrance",
    "observation_error",
    "female_probability",
    "baseline_target_probability",
    "description",
}
MAX_JSON_SAFE_INTEGER = 2**53 - 1
RNG_MODULUS = 1 << 53
RNG_MASK = RNG_MODULUS - 1
RNG_MULTIPLIER = 6_364_136_223_846_793_005
RNG_INCREMENT = 1_442_695_040_888_963_407
MAX_BOOTSTRAP_DRAWS = 5_000_000
MAX_BOOTSTRAP_REFITS = 250_000
MODEL_SNAPSHOT_FIELDS = (
    "mode",
    "normal_allele",
    "variant_allele",
    "normal_label",
    "target_label",
    "intermediate_label",
    "penetrance",
    "observation_error",
    "female_probability",
    "baseline_target_probability",
)


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
    if value != value.strip():
        raise ValidationError(f"{name}_whitespace_forbidden")
    return value


def _integer(
    value: Any,
    name: str,
    *,
    minimum: int,
    maximum: int,
) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not minimum <= value <= maximum
    ):
        raise ValidationError(f"{name}_invalid")
    return value


def _reject_unknown_fields(
    value: Mapping[str, Any],
    allowed: set[str],
    name: str,
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValidationError(f"{name}_unknown_field:{unknown[0]}")


def _derived_seed(base_seed: int, *parts: object) -> int:
    material = "|".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.sha256(material.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big", signed=False) & RNG_MASK


class _CrossRuntimeRandom:
    def __init__(self, seed: int) -> None:
        self._state = seed & RNG_MASK

    def random(self) -> float:
        self._state = (
            self._state * RNG_MULTIPLIER + RNG_INCREMENT
        ) & RNG_MASK
        return self._state / RNG_MODULUS


def _multinomial_counts(
    probabilities: Mapping[str, float],
    total: int,
    rng: _CrossRuntimeRandom,
) -> dict[str, int]:
    categories = sorted(probabilities, key=lambda value: value.encode("utf-8"))
    weights = [float(probabilities[category]) for category in categories]
    probability_total = sum(weights)
    if probability_total <= 0.0:
        raise ValidationError("simulation_probability_total_zero")
    cumulative: list[float] = []
    running = 0.0
    for weight in weights:
        if weight < 0.0 or not math.isfinite(weight):
            raise ValidationError("simulation_probability_invalid")
        running += weight / probability_total
        cumulative.append(running)
    cumulative[-1] = 1.0

    counts = {category: 0 for category in categories}
    for _ in range(total):
        draw = rng.random()
        for index, threshold in enumerate(cumulative):
            if draw < threshold:
                counts[categories[index]] += 1
                break
    return counts


def _apply_observation_error(
    probabilities: Mapping[str, float],
    error_rate: float,
) -> dict[str, float]:
    categories = sorted(probabilities)
    if not categories:
        raise ValidationError("simulation_categories_empty")
    smoothed = {
        category: (
            (1.0 - error_rate) * float(probabilities[category])
            + error_rate / len(categories)
        )
        for category in categories
    }
    total = sum(smoothed.values())
    if total <= 0.0:
        raise ValidationError("simulation_probability_total_zero")
    return {
        category: probability / total
        for category, probability in smoothed.items()
    }


def _observations(counts: Mapping[str, int]) -> list[OffspringObservation]:
    observations: list[OffspringObservation] = []
    for key, count in sorted(counts.items()):
        if count <= 0:
            continue
        try:
            raw_sex, phenotype = key.split(":", 1)
            sex = Sex(raw_sex)
        except (ValueError, KeyError) as exc:
            raise ValidationError(
                f"simulation_phenotype_probability_key_invalid:{key}"
            ) from exc
        observations.append(
            OffspringObservation(
                phenotype=phenotype,
                count=count,
                sex=sex,
            )
        )
    if not observations:
        raise ValidationError("simulation_generated_no_observations")
    return observations


def _order_independent_mean(values: Sequence[float]) -> float:
    # Deterministically sorted Neumaier sum: the mean does not depend on the
    # caller's iteration order and is bit-identical to npm/src/index.mjs
    # orderIndependentMean — keep both sides in sync.
    return _neumaier_sum(sorted(values)) / len(values)


def _entropy_bits(probabilities: Sequence[float]) -> float:
    # Mirrors npm/src/index.mjs entropy (same log2 primitive, same sorted
    # Neumaier accumulation) — keep both sides in sync.
    terms = sorted(
        probability * _portable_log2(probability)
        for probability in probabilities
        if probability > 0.0
    )
    return -_neumaier_sum(terms)


def _quantile(values: Sequence[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _wilson_interval(successes: int, total: int) -> dict[str, float] | None:
    if total <= 0:
        return None
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    radius = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / total
            + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return {
        "low": max(0.0, center - radius),
        "high": min(1.0, center + radius),
    }


def _posterior_discrimination_bins(
    true_posteriors: Sequence[float],
    recovered: Sequence[bool],
) -> list[dict[str, float | int]]:
    bins: list[list[tuple[float, bool]]] = [[] for _ in range(5)]
    for posterior, was_recovered in zip(true_posteriors, recovered):
        index = min(4, int(posterior * 5.0))
        bins[index].append((posterior, was_recovered))
    result: list[dict[str, float | int]] = []
    for index, items in enumerate(bins):
        if not items:
            continue
        result.append(
            {
                "lower": index / 5.0,
                "upper": (index + 1) / 5.0,
                "count": len(items),
                # _order_independent_mean keeps this bit-identical with Node
                # regardless of Python version (builtins.sum gained Neumaier
                # compensation for floats in CPython 3.12).
                "mean_predicted": _order_independent_mean(
                    [item[0] for item in items]
                ),
                "observed_top1_rate": (
                    sum(item[1] for item in items) / len(items)
                ),
            }
        )
    return result


def _deviance(
    counts: Mapping[str, int],
    probabilities: Mapping[str, float],
) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    result = 0.0
    for category, count in counts.items():
        if count <= 0:
            continue
        probability = float(probabilities.get(category, 0.0))
        if probability <= 0.0:
            return float("inf")
        expected = total * probability
        result += 2.0 * count * _portable_log(count / expected)
    return max(0.0, result)


def _goodness_of_fit(
    *,
    models: Sequence[InheritanceModel],
    male_parent_phenotype: str,
    female_parent_phenotype: str,
    pooled_counts: Mapping[str, int],
    bootstrap_replicates: int,
    alpha: float,
    base_seed: int,
    experiment_id: str,
    brood_size: int,
) -> dict[str, Any]:
    observations = _observations(pooled_counts)
    total = sum(pooled_counts.values())
    candidates: list[dict[str, Any]] = []
    for model in models:
        score = score_model(
            model,
            male_parent_phenotype,
            female_parent_phenotype,
            observations,
        )
        configuration = score["best_parent_configuration"]
        if configuration is None:
            candidates.append(
                {
                    "model_id": model.model_id,
                    "deviance": None,
                    "bootstrap_p_value": 0.0,
                    "status": "parent_phenotype_incompatible",
                }
            )
            continue
        predicted = configuration["predicted_probabilities"]
        observed_deviance = _deviance(pooled_counts, predicted)
        seed = _derived_seed(
            base_seed,
            experiment_id,
            brood_size,
            model.model_id,
            "goodness_of_fit",
        )
        rng = _CrossRuntimeRandom(seed)
        extreme = 0
        for _ in range(bootstrap_replicates):
            simulated = _multinomial_counts(predicted, total, rng)
            bootstrap_score = score_model(
                model,
                male_parent_phenotype,
                female_parent_phenotype,
                _observations(simulated),
            )
            bootstrap_configuration = bootstrap_score[
                "best_parent_configuration"
            ]
            bootstrap_deviance = (
                _deviance(
                    simulated,
                    bootstrap_configuration["predicted_probabilities"],
                )
                if bootstrap_configuration is not None
                else float("inf")
            )
            if bootstrap_deviance >= observed_deviance:
                extreme += 1
        p_value = (extreme + 1.0) / (bootstrap_replicates + 1.0)
        candidates.append(
            {
                "model_id": model.model_id,
                "deviance": (
                    observed_deviance
                    if math.isfinite(observed_deviance)
                    else None
                ),
                "bootstrap_p_value": p_value,
                "status": "evaluated",
            }
        )
    all_below = all(
        candidate["bootstrap_p_value"] < alpha for candidate in candidates
    )
    return {
        "method": "parametric_bootstrap_multinomial_deviance",
        "bootstrap_replicates": bootstrap_replicates,
        "alpha": alpha,
        "candidates": candidates,
        "all_candidates_below_alpha": all_below,
    }


def _summarize_replicates(
    *,
    experiment_id: str,
    generating_model_id: str,
    candidate_models: Sequence[InheritanceModel],
    male_parent_phenotype: str,
    female_parent_phenotype: str,
    brood_size: int,
    replicate_results: Sequence[Mapping[str, Any]],
    pooled_counts: Mapping[str, int],
    bootstrap_replicates: int,
    alpha: float,
    base_seed: int,
) -> dict[str, Any]:
    true_in_set = any(
        model.model_id == generating_model_id for model in candidate_models
    )
    recovered = [
        result["best_models"] == [generating_model_id]
        for result in replicate_results
    ]
    top_set_contains_true = [
        generating_model_id in result["best_models"]
        for result in replicate_results
    ]
    true_posteriors = [
        float(result["true_model_posterior"])
        for result in replicate_results
        if result["true_model_posterior"] is not None
    ]
    goodness_of_fit = _goodness_of_fit(
        models=candidate_models,
        male_parent_phenotype=male_parent_phenotype,
        female_parent_phenotype=female_parent_phenotype,
        pooled_counts=pooled_counts,
        bootstrap_replicates=bootstrap_replicates,
        alpha=alpha,
        base_seed=base_seed,
        experiment_id=experiment_id,
        brood_size=brood_size,
    )
    warnings: list[str] = []
    small_sample_warning = brood_size < 30
    if small_sample_warning:
        warnings.append("small_sample_high_uncertainty")
    if goodness_of_fit["all_candidates_below_alpha"]:
        warnings.append("candidate_model_set_misspecified_suspected")

    best_model_set_counts: dict[str, int] = defaultdict(int)
    best_model_membership_counts: dict[str, int] = defaultdict(int)
    unique_best_model_counts: dict[str, int] = defaultdict(int)
    for result in replicate_results:
        best_models = result["best_models"]
        best_model_set_counts["|".join(best_models)] += 1
        for model_id in best_models:
            best_model_membership_counts[model_id] += 1
        if len(best_models) == 1:
            unique_best_model_counts[best_models[0]] += 1
    replicate_count = len(replicate_results)
    successes = sum(recovered) if true_in_set else 0
    brier_values = [
        float(result["multiclass_brier"])
        for result in replicate_results
        if result["multiclass_brier"] is not None
    ]
    return {
        "brood_size": brood_size,
        "replicates": replicate_count,
        "small_sample_warning": small_sample_warning,
        "warnings": warnings,
        "recovery": {
            "true_model_in_candidate_set": true_in_set,
            "top1_rate": successes / replicate_count if true_in_set else None,
            "top_set_contains_true_rate": (
                sum(top_set_contains_true) / replicate_count
                if true_in_set
                else None
            ),
            "tie_policy": "top1_rate_requires_unique_best_model",
            "wilson_95": (
                _wilson_interval(successes, replicate_count)
                if true_in_set
                else None
            ),
            "mean_true_model_posterior": (
                _order_independent_mean(true_posteriors)
                if true_posteriors
                else None
            ),
            "multiclass_brier_mean": (
                _order_independent_mean(brier_values)
                if brier_values
                else None
            ),
        },
        "uncertainty": {
            "true_model_posterior_quantiles": {
                "p05": _quantile(true_posteriors, 0.05),
                "p50": _quantile(true_posteriors, 0.50),
                "p95": _quantile(true_posteriors, 0.95),
            },
            "mean_model_entropy_bits": _order_independent_mean(
                [
                    float(result["model_entropy_bits"])
                    for result in replicate_results
                ]
            ),
            "top_two_log2_bayes_factor_quantiles": {
                "p05": _quantile(
                    [
                        float(result["top_two_log2_bayes_factor"])
                        for result in replicate_results
                        if result["top_two_log2_bayes_factor"] is not None
                    ],
                    0.05,
                ),
                "p50": _quantile(
                    [
                        float(result["top_two_log2_bayes_factor"])
                        for result in replicate_results
                        if result["top_two_log2_bayes_factor"] is not None
                    ],
                    0.50,
                ),
                "p95": _quantile(
                    [
                        float(result["top_two_log2_bayes_factor"])
                        for result in replicate_results
                        if result["top_two_log2_bayes_factor"] is not None
                    ],
                    0.95,
                ),
            },
        },
        "best_model_frequencies": {
            model_id: {
                "top_set_count": best_model_membership_counts[model_id],
                "top_set_rate": (
                    best_model_membership_counts[model_id] / replicate_count
                ),
                "unique_top1_count": unique_best_model_counts.get(
                    model_id,
                    0,
                ),
                "unique_top1_rate": (
                    unique_best_model_counts.get(model_id, 0)
                    / replicate_count
                ),
            }
            for model_id in sorted(best_model_membership_counts)
        },
        "best_model_set_frequencies": {
            model_set: {
                "count": count,
                "rate": count / replicate_count,
            }
            for model_set, count in sorted(best_model_set_counts.items())
        },
        "posterior_discrimination_bins": (
            _posterior_discrimination_bins(true_posteriors, recovered)
            if true_in_set
            else []
        ),
        "goodness_of_fit": goodness_of_fit,
        "replicate_results": list(replicate_results),
    }


def _parse_experiments(
    value: Any,
) -> list[dict[str, Any]]:
    raw_experiments = _sequence(value, "simulation_experiments")
    if not raw_experiments:
        raise ValidationError("simulation_experiments_required")
    available_models = {
        model.model_id: model for model in default_candidate_models()
    }
    parsed: list[dict[str, Any]] = []
    seen_experiment_ids: set[str] = set()
    total_inference_runs = 0
    total_simulated_offspring = 0
    for raw in raw_experiments:
        item = _mapping(raw, "simulation_experiment")
        _reject_unknown_fields(item, EXPERIMENT_FIELDS, "simulation_experiment")
        experiment_id = _required_text(
            item.get("experiment_id"),
            "simulation_experiment_id",
        )
        if experiment_id in seen_experiment_ids:
            raise ValidationError(
                f"simulation_experiment_id_duplicate:{experiment_id}"
            )
        seen_experiment_ids.add(experiment_id)

        raw_model = _mapping(
            item.get("generating_model"),
            "simulation_generating_model",
        )
        _reject_unknown_fields(
            raw_model,
            MODEL_FIELDS,
            "simulation_generating_model",
        )
        _required_text(
            raw_model.get("model_id"),
            "simulation_generating_model_id",
        )
        _required_text(
            raw_model.get("mode"),
            "simulation_generating_model_mode",
        )
        for field in (
            "normal_allele",
            "variant_allele",
            "normal_label",
            "target_label",
            "intermediate_label",
        ):
            if field in raw_model:
                _required_text(
                    raw_model[field],
                    f"simulation_generating_model_{field}",
                )
        if "description" in raw_model:
            description = raw_model["description"]
            if not isinstance(description, str):
                raise ValidationError(
                    "simulation_generating_model_description_invalid"
                )
            if description != description.strip():
                raise ValidationError(
                    "simulation_generating_model_description_"
                    "whitespace_forbidden"
                )
        for field in (
            "penetrance",
            "observation_error",
            "female_probability",
            "baseline_target_probability",
        ):
            if field not in raw_model:
                continue
            probability = raw_model[field]
            if (
                isinstance(probability, bool)
                or not isinstance(probability, (int, float))
                or not math.isfinite(float(probability))
                or not 0.0 <= float(probability) <= 1.0
            ):
                raise ValidationError(
                    f"simulation_generating_model_{field}_invalid"
                )
        generating_model = InheritanceModel.from_dict(raw_model)
        raw_candidate_ids = _sequence(
            item.get("candidate_model_ids"),
            "simulation_candidate_model_ids",
        )
        candidate_ids = [
            _required_text(value, "simulation_candidate_model_id")
            for value in raw_candidate_ids
        ]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValidationError(
                f"simulation_candidate_model_id_duplicate:{experiment_id}"
            )
        if candidate_ids.count("M0_NONGENETIC") != 1:
            raise ValidationError(
                "simulation_candidate_set_requires_exactly_one_"
                f"M0_NONGENETIC:{experiment_id}"
            )
        if len(candidate_ids) < 2:
            raise ValidationError(
                "simulation_candidate_set_requires_multiple_models:"
                f"{experiment_id}"
            )
        unknown_models = sorted(set(candidate_ids) - set(available_models))
        if unknown_models:
            raise ValidationError(
                f"simulation_candidate_model_unknown:{unknown_models[0]}"
            )
        candidate_models = [
            available_models[model_id]
            for model_id in candidate_ids
        ]
        candidate_models = [
            InheritanceModel(
                **{
                    **model.to_dict(),
                    "mode": model.mode,
                    "target_label": generating_model.target_label,
                    "normal_label": generating_model.normal_label,
                    "intermediate_label": generating_model.intermediate_label,
                }
            )
            for model in candidate_models
        ]
        if generating_model.model_id in available_models:
            canonical_model = available_models[generating_model.model_id]
            matching_candidate = InheritanceModel(
                **{
                    **canonical_model.to_dict(),
                    "mode": canonical_model.mode,
                    "target_label": generating_model.target_label,
                    "normal_label": generating_model.normal_label,
                    "intermediate_label": generating_model.intermediate_label,
                }
            )
            if any(
                getattr(generating_model, field)
                != getattr(matching_candidate, field)
                for field in MODEL_SNAPSHOT_FIELDS
            ):
                raise ValidationError(
                    "simulation_generating_model_snapshot_mismatch:"
                    f"{experiment_id}:{generating_model.model_id}"
                )

        raw_brood_sizes = _sequence(
            item.get("brood_sizes"),
            "simulation_brood_sizes",
        )
        brood_sizes = [
            _integer(
                brood_size,
                "simulation_brood_size",
                minimum=1,
                maximum=100_000,
            )
            for brood_size in raw_brood_sizes
        ]
        if not brood_sizes:
            raise ValidationError(
                f"simulation_brood_sizes_required:{experiment_id}"
            )
        if len(brood_sizes) != len(set(brood_sizes)):
            raise ValidationError(
                f"simulation_brood_size_duplicate:{experiment_id}"
            )
        brood_sizes.sort()
        replicates = _integer(
            item.get("replicates"),
            "simulation_replicates",
            minimum=1,
            maximum=10_000,
        )
        male_genotype = _required_text(
            item.get("male_genotype"),
            "simulation_male_genotype",
        )
        female_genotype = _required_text(
            item.get("female_genotype"),
            "simulation_female_genotype",
        )
        male_parent_phenotype = _required_text(
            item.get("male_parent_phenotype"),
            "simulation_male_parent_phenotype",
        )
        female_parent_phenotype = _required_text(
            item.get("female_parent_phenotype"),
            "simulation_female_parent_phenotype",
        )
        if generating_model.mode.value == "non_genetic":
            for role, genotype in (
                ("male", male_genotype),
                ("female", female_genotype),
            ):
                if genotype != "NA":
                    raise ValidationError(
                        "simulation_non_genetic_genotype_must_be_NA:"
                        f"{experiment_id}:{role}"
                    )
        for role, sex, genotype, parent_phenotype in (
            (
                "male",
                Sex.MALE,
                male_genotype,
                male_parent_phenotype,
            ),
            (
                "female",
                Sex.FEMALE,
                female_genotype,
                female_parent_phenotype,
            ),
        ):
            parent_distribution = phenotype_distribution(
                generating_model,
                sex,
                genotype,
            )
            if (
                parent_phenotype not in {"unknown", "UNK"}
                and parent_distribution.get(parent_phenotype, 0.0)
                <= 0.0
            ):
                raise ValidationError(
                    "simulation_parent_phenotype_genotype_incompatible:"
                    f"{experiment_id}:{role}"
                )
        total_inference_runs += len(brood_sizes) * replicates
        total_simulated_offspring += sum(brood_sizes) * replicates
        parsed.append(
            {
                "experiment_id": experiment_id,
                "generating_model": generating_model,
                "male_genotype": male_genotype,
                "female_genotype": female_genotype,
                "male_parent_phenotype": male_parent_phenotype,
                "female_parent_phenotype": female_parent_phenotype,
                "brood_sizes": brood_sizes,
                "replicates": replicates,
                "candidate_model_ids": candidate_ids,
                "candidate_models": candidate_models,
            }
        )
    if total_inference_runs > 10_000 or total_simulated_offspring > 5_000_000:
        raise ValidationError("simulation_work_budget_exceeded")
    return parsed


def _reliability_bins(
    values: Sequence[tuple[float, bool]],
) -> list[dict[str, float | int]]:
    bins: list[list[tuple[float, bool]]] = [[] for _ in range(5)]
    for posterior, is_true_model in values:
        bins[min(4, int(posterior * 5.0))].append(
            (posterior, is_true_model)
        )
    result: list[dict[str, float | int]] = []
    for index, items in enumerate(bins):
        if not items:
            continue
        result.append(
            {
                "lower": index / 5.0,
                "upper": (index + 1) / 5.0,
                "count": len(items),
                "mean_posterior": _order_independent_mean(
                    [item[0] for item in items]
                ),
                "observed_true_model_rate": (
                    sum(item[1] for item in items) / len(items)
                ),
            }
        )
    return result


def _calibration_by_design(
    experiments: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    groups: dict[
        tuple[int, tuple[str, ...]],
        list[tuple[str, Mapping[str, float]]],
    ] = defaultdict(list)
    for experiment in experiments:
        true_model_id = str(experiment["generating_model_id"])
        candidate_ids = tuple(experiment["candidate_model_ids"])
        for summary in experiment["by_brood_size"]:
            key = (int(summary["brood_size"]), candidate_ids)
            for replicate in summary["replicate_results"]:
                groups[key].append(
                    (
                        true_model_id,
                        replicate["model_posteriors"],
                    )
                )

    result: list[dict[str, Any]] = []
    for (brood_size, candidate_ids), records in groups.items():
        true_model_counts: dict[str, int] = defaultdict(int)
        for true_model_id, _posterior in records:
            true_model_counts[true_model_id] += 1
        true_models = set(true_model_counts)
        status = "estimated"
        reason = None
        if len(true_models) < 2:
            status = "not_estimable"
            reason = "requires_at_least_two_generating_models"
        elif not true_models.issubset(set(candidate_ids)):
            status = "not_estimable"
            reason = "generating_model_outside_candidate_set"

        models: list[dict[str, Any]] = []
        if status == "estimated":
            for model_id in candidate_ids:
                values = [
                    (
                        float(posterior[model_id]),
                        true_model_id == model_id,
                    )
                    for true_model_id, posterior in records
                ]
                models.append(
                    {
                        "model_id": model_id,
                        "bins": _reliability_bins(values),
                    }
                )
        result.append(
            {
                "brood_size": brood_size,
                "candidate_model_ids": list(candidate_ids),
                "method": (
                    "empirical_reliability_under_declared_simulation_design"
                ),
                "status": status,
                "reason": reason,
                "design_model_frequencies": dict(
                    sorted(true_model_counts.items())
                ),
                "models": models,
                "limits": {
                    "design_dependent": True,
                    "not_external_field_calibration": True,
                },
            }
        )
    return result


def simulate_recovery(request: Mapping[str, Any]) -> dict[str, Any]:
    request = _mapping(request, "simulation_request")
    _reject_unknown_fields(request, ROOT_FIELDS, "simulation_request")
    schema_version = request.get("schema_version")
    if schema_version != "1.0.0":
        raise ValidationError(
            f"simulation_schema_version_unsupported:{schema_version}"
        )
    base_seed = _integer(
        request.get("base_seed"),
        "simulation_base_seed",
        minimum=0,
        maximum=MAX_JSON_SAFE_INTEGER,
    )
    bootstrap_replicates = _integer(
        request.get("goodness_of_fit_replicates"),
        "goodness_of_fit_replicates",
        minimum=10,
        maximum=10_000,
    )
    alpha_value = request.get("misspecification_alpha")
    if (
        isinstance(alpha_value, bool)
        or not isinstance(alpha_value, (int, float))
        or not math.isfinite(float(alpha_value))
        or not 0.0 < float(alpha_value) <= 0.25
    ):
        raise ValidationError("misspecification_alpha_invalid")
    alpha = float(alpha_value)
    if bootstrap_replicates + 1 <= 1.0 / alpha:
        raise ValidationError(
            "goodness_of_fit_replicates_insufficient_for_alpha:"
            f"{bootstrap_replicates}:{alpha}"
        )
    experiments = _parse_experiments(request.get("experiments"))
    bootstrap_draws = sum(
        sum(experiment["brood_sizes"])
        * experiment["replicates"]
        * bootstrap_replicates
        * len(experiment["candidate_models"])
        for experiment in experiments
    )
    bootstrap_refits = sum(
        len(experiment["brood_sizes"])
        * bootstrap_replicates
        * len(experiment["candidate_models"])
        for experiment in experiments
    )
    if (
        bootstrap_draws > MAX_BOOTSTRAP_DRAWS
        or bootstrap_refits > MAX_BOOTSTRAP_REFITS
    ):
        raise ValidationError("simulation_bootstrap_work_budget_exceeded")

    output_experiments: list[dict[str, Any]] = []
    for experiment in experiments:
        generating_model = experiment["generating_model"]
        cross_result = cross_known(
            generating_model,
            experiment["male_genotype"],
            experiment["female_genotype"],
        )
        generating_probabilities = _apply_observation_error(
            cross_result["phenotype_probabilities"],
            generating_model.observation_error,
        )
        by_brood_size: list[dict[str, Any]] = []
        for brood_size in experiment["brood_sizes"]:
            replicate_results: list[dict[str, Any]] = []
            pooled_counts: dict[str, int] = defaultdict(int)
            for replicate_index in range(experiment["replicates"]):
                seed = _derived_seed(
                    base_seed,
                    experiment["experiment_id"],
                    brood_size,
                    replicate_index,
                )
                counts = _multinomial_counts(
                    generating_probabilities,
                    brood_size,
                    _CrossRuntimeRandom(seed),
                )
                for category, count in counts.items():
                    pooled_counts[category] += count
                inference = infer_candidate_models(
                    experiment["candidate_models"],
                    experiment["male_parent_phenotype"],
                    experiment["female_parent_phenotype"],
                    _observations(counts),
                )
                posterior = {
                    score["model_id"]: float(score["posterior"])
                    for score in inference["ranking"]
                }
                ranking = inference["ranking"]
                maximum_posterior = max(posterior.values())
                best_models = sorted(
                    model_id
                    for model_id, probability in posterior.items()
                    if math.isclose(
                        probability,
                        maximum_posterior,
                        rel_tol=0.0,
                        abs_tol=1e-12,
                    )
                )
                unique_best_model = (
                    best_models[0] if len(best_models) == 1 else None
                )
                top_two_log2_bayes_factor = None
                if (
                    len(ranking) >= 2
                    and math.isfinite(float(ranking[0]["log_evidence"]))
                    and math.isfinite(float(ranking[1]["log_evidence"]))
                ):
                    top_two_log2_bayes_factor = (
                        float(ranking[0]["log_evidence"])
                        - float(ranking[1]["log_evidence"])
                    ) / _LN2
                true_model_id = generating_model.model_id
                true_posterior = posterior.get(true_model_id)
                brier = None
                if true_posterior is not None:
                    # Sorted Neumaier sum keeps this bit-identical with the
                    # Node implementation on every platform/Python version.
                    brier = _neumaier_sum(
                        sorted(
                            (
                                probability
                                - (1.0 if model_id == true_model_id else 0.0)
                            )
                            ** 2
                            for model_id, probability in posterior.items()
                        )
                    )
                replicate_results.append(
                    {
                        "replicate_index": replicate_index,
                        "seed": seed,
                        "observed_counts": dict(sorted(counts.items())),
                        "best_model": unique_best_model,
                        "best_models": best_models,
                        "top1_tied": len(best_models) > 1,
                        "best_model_posterior": maximum_posterior,
                        "model_posteriors": dict(sorted(posterior.items())),
                        "true_model_posterior": true_posterior,
                        "model_entropy_bits": _entropy_bits(
                            list(posterior.values())
                        ),
                        "top_two_log2_bayes_factor": (
                            top_two_log2_bayes_factor
                        ),
                        "multiclass_brier": brier,
                    }
                )
            by_brood_size.append(
                _summarize_replicates(
                    experiment_id=experiment["experiment_id"],
                    generating_model_id=generating_model.model_id,
                    candidate_models=experiment["candidate_models"],
                    male_parent_phenotype=experiment[
                        "male_parent_phenotype"
                    ],
                    female_parent_phenotype=experiment[
                        "female_parent_phenotype"
                    ],
                    brood_size=brood_size,
                    replicate_results=replicate_results,
                    pooled_counts=dict(pooled_counts),
                    bootstrap_replicates=bootstrap_replicates,
                    alpha=alpha,
                    base_seed=base_seed,
                )
            )
        output_experiments.append(
            {
                "experiment_id": experiment["experiment_id"],
                "generating_model_id": generating_model.model_id,
                "candidate_model_ids": experiment["candidate_model_ids"],
                "by_brood_size": by_brood_size,
            }
        )
    return {
        "schema_version": schema_version,
        "base_seed": base_seed,
        "experiments": output_experiments,
        "calibration_by_design": _calibration_by_design(
            output_experiments
        ),
        "evidence": {
            "source": "synthetic_simulation_only",
            "evidence_grade": "G0",
            "gene_status": "not_assessed",
        },
        "limits": {
            "synthetic_only": True,
            "decision_support_only": True,
            "posterior_scope": "relative_to_candidate_model_set",
            "gene_status": "not_assessed",
            "inheritance_status": "simulation_recovery_not_confirmation",
            "genotype_specific_survival": "not_modeled",
        },
    }
