from __future__ import annotations

import math
import struct
from collections import defaultdict
from itertools import product
from typing import Any, Iterable, Mapping, Sequence

from .cross import collapse_phenotypes_by_sex, cross_known
from .enums import InheritanceMode, Sex
from .errors import ValidationError
from .genotype import compatible_genotypes
from .models import InheritanceModel, OffspringObservation


# Cross-runtime numeric primitives. Each helper below is mirrored
# operation-for-operation in npm/src/index.mjs (portableLog, portableExp,
# neumaierSum, logSumExp, logGamma) so both runtimes produce bit-identical
# floats by construction — keep the two implementations in sync.
#
# Why not the platform libm: glibc and V8 disagree in the last bit for
# roughly 3-10% of log/exp inputs (measured on Linux x86-64), so any
# likelihood built on math.log/Math.log can never be exactly reproducible
# across runtimes. The fdlibm algorithms below use only IEEE-754
# correctly-rounded operations (+, -, *, /) plus explicit bit manipulation,
# which both runtimes define identically.
_DOUBLE = struct.Struct("<d")
_UINT64 = struct.Struct("<Q")


def _double_words(value: float) -> tuple[int, int]:
    bits = _UINT64.unpack(_DOUBLE.pack(value))[0]
    return bits >> 32, bits & 0xFFFFFFFF


def _words_double(high: int, low: int) -> float:
    return _DOUBLE.unpack(
        _UINT64.pack(((high & 0xFFFFFFFF) << 32) | (low & 0xFFFFFFFF))
    )[0]


_LOG_LN2_HI = 6.93147180369123816490e-01
_LOG_LN2_LO = 1.90821492927058770002e-10
_LOG_LG1 = 6.666666666666735130e-01
_LOG_LG2 = 3.999999999940941908e-01
_LOG_LG3 = 2.857142874366239149e-01
_LOG_LG4 = 2.222219843214978396e-01
_LOG_LG5 = 1.818357216161805012e-01
_LOG_LG6 = 1.531383769920937332e-01
_LOG_LG7 = 1.479819860511658591e-01
_TWO54 = 1.8014398509481984e16


def _portable_log(x: float) -> float:
    # fdlibm e_log.c (Sun Microsystems), ported verbatim; mirrored in
    # npm/src/index.mjs portableLog.
    high, low = _double_words(x)
    k = 0
    if high < 0x00100000 or high >= 0x80000000:
        if ((high & 0x7FFFFFFF) | low) == 0:
            return float("-inf")
        if high >= 0x80000000:
            return float("nan")
        k -= 54
        x *= _TWO54
        high, low = _double_words(x)
    if high >= 0x7FF00000:
        return x + x
    k += (high >> 20) - 1023
    high &= 0x000FFFFF
    i = (high + 0x95F64) & 0x100000
    x = _words_double(high | (i ^ 0x3FF00000), low)
    k += i >> 20
    f = x - 1.0
    if (0x000FFFFF & (2 + high)) < 3:
        if f == 0.0:
            if k == 0:
                return 0.0
            dk = float(k)
            return dk * _LOG_LN2_HI + dk * _LOG_LN2_LO
        r = f * f * (0.5 - 0.33333333333333333 * f)
        if k == 0:
            return f - r
        dk = float(k)
        return dk * _LOG_LN2_HI - ((r - dk * _LOG_LN2_LO) - f)
    s = f / (2.0 + f)
    dk = float(k)
    z = s * s
    i = high - 0x6147A
    w = z * z
    j = 0x6B851 - high
    t1 = w * (_LOG_LG2 + w * (_LOG_LG4 + w * _LOG_LG6))
    t2 = z * (_LOG_LG1 + w * (_LOG_LG3 + w * (_LOG_LG5 + w * _LOG_LG7)))
    i |= j
    r = t2 + t1
    if i > 0:
        hfsq = 0.5 * f * f
        if k == 0:
            return f - (hfsq - s * (hfsq + r))
        return dk * _LOG_LN2_HI - (
            (hfsq - (s * (hfsq + r) + dk * _LOG_LN2_LO)) - f
        )
    if k == 0:
        return f - s * (f - r)
    return dk * _LOG_LN2_HI - ((s * (f - r) - dk * _LOG_LN2_LO) - f)


_EXP_P1 = 1.66666666666666019037e-01
_EXP_P2 = -2.77777777770155933842e-03
_EXP_P3 = 6.61375632143793436117e-05
_EXP_P4 = -1.65339022054652515390e-06
_EXP_P5 = 4.13813679705723846039e-08
_EXP_LN2_HI = 6.93147180369123816490e-01
_EXP_LN2_LO = 1.90821492927058770002e-10
_EXP_INV_LN2 = 1.44269504088896338700e00
_EXP_O_THRESHOLD = 7.09782712893383973096e02
_EXP_U_THRESHOLD = -7.45133219101941108420e02
_EXP_HUGE = 1.0e300
_EXP_TWOM1000 = 9.33263618503218878990e-302


def _portable_exp(x: float) -> float:
    # fdlibm e_exp.c (Sun Microsystems), ported verbatim; mirrored in
    # npm/src/index.mjs portableExp.
    high, low = _double_words(x)
    xsb = (high >> 31) & 1
    hx = high & 0x7FFFFFFF
    if hx >= 0x40862E42:
        if hx >= 0x7FF00000:
            if ((hx & 0xFFFFF) | low) != 0:
                return x + x
            return 0.0 if xsb == 1 else x
        if x > _EXP_O_THRESHOLD:
            return _EXP_HUGE * _EXP_HUGE
        if x < _EXP_U_THRESHOLD:
            return _EXP_TWOM1000 * _EXP_TWOM1000
    k = 0
    hi = 0.0
    lo = 0.0
    if hx > 0x3FD62E42:
        if hx < 0x3FF0A2E4:
            hi = x - (_EXP_LN2_HI if xsb == 0 else -_EXP_LN2_HI)
            lo = _EXP_LN2_LO if xsb == 0 else -_EXP_LN2_LO
            k = 1 - xsb - xsb
        else:
            k = int(_EXP_INV_LN2 * x + (0.5 if xsb == 0 else -0.5))
            t = float(k)
            hi = x - t * _EXP_LN2_HI
            lo = t * _EXP_LN2_LO
        x = hi - lo
    elif hx < 0x3E300000:
        if _EXP_HUGE + x > 1.0:
            return 1.0 + x
    t = x * x
    c = x - t * (
        _EXP_P1
        + t * (_EXP_P2 + t * (_EXP_P3 + t * (_EXP_P4 + t * _EXP_P5)))
    )
    if k == 0:
        return 1.0 - ((x * c / (c - 2.0)) - x)
    y = 1.0 - ((lo - (x * c) / (2.0 - c)) - hi)
    high_y, low_y = _double_words(y)
    if k >= -1021:
        return _words_double(high_y + (k << 20), low_y)
    return _words_double(high_y + ((k + 1000) << 20), low_y) * _EXP_TWOM1000


_LN2 = _portable_log(2.0)
_INV_LN2 = 1.4426950408889634


def _portable_log2(value: float) -> float:
    # Deterministic (not correctly-rounded) log2: same expression in both
    # runtimes, so bit-identical by construction.
    return _portable_log(value) * _INV_LN2
def _neumaier_sum(values: Sequence[float]) -> float:
    # Neumaier-compensated summation. Only IEEE-exact operations (+, -, abs,
    # comparison), so the result is bit-identical across runtimes for the
    # same input order. Callers sort summands first to also make the result
    # independent of the caller's iteration order.
    total = 0.0
    compensation = 0.0
    for value in values:
        candidate = total + value
        if abs(total) >= abs(value):
            compensation += (total - candidate) + value
        else:
            compensation += (value - candidate) + total
        total = candidate
    return total + compensation


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        return float("-inf")
    maximum = max(values)
    if math.isinf(maximum):
        return maximum
    # Sorting the exponentiated terms makes the sum independent of the
    # caller's model/configuration ordering (order-independence guarantee).
    terms = sorted(_portable_exp(value - maximum) for value in values)
    return maximum + _portable_log(_neumaier_sum(terms))


_LANCZOS_COEFFICIENTS = (
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7,
)


def _log_gamma(value: float) -> float:
    # Lanczos approximation shared verbatim with npm/src/index.mjs logGamma
    # so Python and Node run one algorithm (instead of math.lgamma vs a
    # hand-rolled series, which agree only to ~1e-13 relative). Built on the
    # portable log so it is bit-identical across runtimes for value >= 0.5;
    # engine call sites only evaluate integers >= 1. The reflection branch
    # (value < 0.5) keeps native math.sin and is not bit-guaranteed.
    if value < 0.5:
        return (
            _portable_log(math.pi)
            - _portable_log(math.sin(math.pi * value))
            - _log_gamma(1.0 - value)
        )
    shifted = value - 1.0
    series = _LANCZOS_COEFFICIENTS[0]
    for index in range(1, len(_LANCZOS_COEFFICIENTS)):
        series += _LANCZOS_COEFFICIENTS[index] / (shifted + index)
    scale = shifted + 7.5
    return (
        0.5 * _portable_log(2.0 * math.pi)
        + (shifted + 0.5) * _portable_log(scale)
        - scale
        + _portable_log(series)
    )


def _normalize_log_weights(values: Mapping[str, float]) -> dict[str, float]:
    denominator = _logsumexp(list(values.values()))
    if math.isinf(denominator):
        uniform = 1.0 / len(values)
        return {key: uniform for key in values}
    return {
        key: _portable_exp(value - denominator)
        for key, value in values.items()
    }


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
    coefficient = _log_gamma(total + 1.0) - sum(
        _log_gamma(count + 1.0) for count in counts.values()
    )
    likelihood = coefficient
    for category, count in counts.items():
        probability = float(probabilities.get(category, 0.0))
        if count > 0 and probability <= 0.0:
            return float("-inf")
        if count > 0:
            likelihood += count * _portable_log(probability)
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
        log_joint = _portable_log(prior) + log_likelihood
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
            _portable_exp(configuration["log_joint"] - log_evidence)
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


CONFIGURATION_PRIOR_WARNING = (
    "model_ranking_may_reflect_uniform_configuration_prior_not_data"
)


def _configuration_prior_warnings(
    scores: Sequence[Mapping[str, Any]],
) -> list[str]:
    # Statistical caution, not a biological conclusion. score_model() gives
    # every parent configuration of a model the uniform prior 1/n (n =
    # len(configurations)), so its log evidence decomposes exactly as
    #     log_evidence = logsumexp_c(log_likelihood_c) - log(n).
    # For the top-ranked model T and any runner-up R the evidence gap is
    #     log_evidence_T - log_evidence_R
    #         = [logsumexp(ll_T) - logsumexp(ll_R)] + [log(n_R) - log(n_T)]
    # where the bracketed second term is contributed by configuration counts
    # (model complexity) alone. Whenever
    #     log_evidence_T - log_evidence_R <= log(n_R) - log(n_T),
    # equivalently
    #     log_evidence_T + log(n_T) <= log_evidence_R + log(n_R),
    # the flat (prior-free) configuration likelihood of R is at least that of
    # T, so the ranking gap could be explained by the uniform configuration
    # prior rather than the offspring data. Mirrored in npm/src/index.mjs
    # (inferCandidateModels) — keep both sides in sync.
    top = scores[0]
    top_configurations = len(top["configurations"])
    if top_configurations == 0 or not math.isfinite(top["log_evidence"]):
        return []
    top_flat = top["log_evidence"] + _portable_log(top_configurations)
    for other in scores[1:]:
        other_configurations = len(other["configurations"])
        if other_configurations == 0 or not math.isfinite(
            other["log_evidence"]
        ):
            continue
        if top_flat <= other["log_evidence"] + _portable_log(
            other_configurations
        ):
            return [CONFIGURATION_PRIOR_WARNING]
    return []


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
    # Z-linked candidates always predict a 0.5 offspring sex ratio (maternal
    # ZW segregation), so mixing them with candidates that assume a different
    # female_probability would bias model comparison silently. Fail closed
    # (mirrored in npm/src/index.mjs inferCandidateModels — keep in sync).
    if any(
        model.mode
        in {
            InheritanceMode.Z_LINKED_RECESSIVE,
            InheritanceMode.Z_LINKED_DOMINANT,
        }
        for model in models
    ) and any(model.female_probability != 0.5 for model in models):
        raise ValidationError(
            "female_probability_must_be_0_5_when_z_linked_model_present"
        )

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
            score["log_evidence"] + _portable_log(priors[model.model_id])
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
        "warnings": _configuration_prior_warnings(scores),
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
