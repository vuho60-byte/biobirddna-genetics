from __future__ import annotations

import math
from collections import Counter
from typing import Iterable

from .enums import GenotypeState, InheritanceMode, Sex
from .errors import UnsupportedModelError, ValidationError
from .models import InheritanceModel


def _canonical_pair(a: str, b: str, model: InheritanceModel) -> str:
    order = {model.normal_allele: 0, model.variant_allele: 1}
    try:
        left, right = sorted((a, b), key=lambda allele: order[allele])
    except KeyError as exc:
        raise ValidationError(f"unknown_allele:{exc.args[0]}") from exc
    return f"{left}/{right}"


def parse_genotype(model: InheritanceModel, sex: Sex, value: str) -> str:
    """Validate and canonicalize the compact genotype representation.

    Autosomal and Z-linked males use two allele tokens, e.g. ``N/m``.
    Z-linked females use one Z allele plus W, e.g. ``N/W`` or ``m/W``.
    """

    if not isinstance(value, str) or not value.strip():
        raise ValidationError("genotype_required")
    tokens = [token.strip() for token in value.split("/")]
    allowed = {model.normal_allele, model.variant_allele}

    if model.mode == InheritanceMode.NON_GENETIC:
        if value != "NA":
            raise ValidationError("non_genetic_model_requires_NA_genotype")
        return "NA"

    if model.mode in {
        InheritanceMode.AUTOSOMAL_RECESSIVE,
        InheritanceMode.AUTOSOMAL_DOMINANT,
        InheritanceMode.INCOMPLETE_DOMINANT,
    }:
        if len(tokens) != 2 or any(token not in allowed for token in tokens):
            raise ValidationError("autosomal_genotype_invalid")
        return _canonical_pair(tokens[0], tokens[1], model)

    if model.mode in {
        InheritanceMode.Z_LINKED_RECESSIVE,
        InheritanceMode.Z_LINKED_DOMINANT,
    }:
        if sex == Sex.MALE:
            if len(tokens) != 2 or any(token not in allowed for token in tokens):
                raise ValidationError("z_linked_male_genotype_invalid")
            return _canonical_pair(tokens[0], tokens[1], model)
        if sex == Sex.FEMALE:
            if len(tokens) != 2 or tokens[1] != "W" or tokens[0] not in allowed:
                raise ValidationError("z_linked_female_genotype_invalid")
            return f"{tokens[0]}/W"
        raise ValidationError("sex_required_for_z_linked_genotype")

    raise UnsupportedModelError(f"unsupported_mode:{model.mode.value}")


def enumerate_genotypes(model: InheritanceModel, sex: Sex) -> list[str]:
    n = model.normal_allele
    m = model.variant_allele

    if model.mode == InheritanceMode.NON_GENETIC:
        return ["NA"]
    if model.mode in {
        InheritanceMode.AUTOSOMAL_RECESSIVE,
        InheritanceMode.AUTOSOMAL_DOMINANT,
        InheritanceMode.INCOMPLETE_DOMINANT,
    }:
        return [f"{n}/{n}", f"{n}/{m}", f"{m}/{m}"]
    if model.mode in {
        InheritanceMode.Z_LINKED_RECESSIVE,
        InheritanceMode.Z_LINKED_DOMINANT,
    }:
        if sex == Sex.MALE:
            return [f"{n}/{n}", f"{n}/{m}", f"{m}/{m}"]
        if sex == Sex.FEMALE:
            return [f"{n}/W", f"{m}/W"]
        raise ValidationError("sex_required_for_z_linked_genotype")
    raise UnsupportedModelError(f"unsupported_mode:{model.mode.value}")


def allele_counts(genotype: str) -> Counter[str]:
    return Counter(token for token in genotype.split("/") if token != "W")


def genotype_state(
    model: InheritanceModel,
    sex: Sex,
    genotype: str,
) -> GenotypeState:
    genotype = parse_genotype(model, sex, genotype)
    if model.mode == InheritanceMode.NON_GENETIC:
        return GenotypeState.UNKNOWN

    copies = allele_counts(genotype)[model.variant_allele]

    if model.mode in {
        InheritanceMode.AUTOSOMAL_RECESSIVE,
        InheritanceMode.Z_LINKED_RECESSIVE,
    }:
        required = 1 if (model.mode == InheritanceMode.Z_LINKED_RECESSIVE and sex == Sex.FEMALE) else 2
        if copies >= required:
            return GenotypeState.AFFECTED
        if copies == 1:
            return GenotypeState.CARRIER
        return GenotypeState.NORMAL

    if model.mode in {
        InheritanceMode.AUTOSOMAL_DOMINANT,
        InheritanceMode.Z_LINKED_DOMINANT,
    }:
        return GenotypeState.AFFECTED if copies >= 1 else GenotypeState.NORMAL

    if model.mode == InheritanceMode.INCOMPLETE_DOMINANT:
        if copies == 0:
            return GenotypeState.NORMAL
        if copies == 1:
            return GenotypeState.INTERMEDIATE
        return GenotypeState.AFFECTED

    raise UnsupportedModelError(f"unsupported_mode:{model.mode.value}")


def phenotype_distribution(
    model: InheritanceModel,
    sex: Sex,
    genotype: str,
) -> dict[str, float]:
    """Return visible phenotype probabilities after penetrance is applied."""

    if model.mode == InheritanceMode.NON_GENETIC:
        p = model.baseline_target_probability
        return {model.target_label: p, model.normal_label: 1.0 - p}

    state = genotype_state(model, sex, genotype)
    if state == GenotypeState.INTERMEDIATE:
        return {model.intermediate_label: 1.0}
    if state == GenotypeState.AFFECTED:
        p = model.penetrance
        return {model.target_label: p, model.normal_label: 1.0 - p}
    return {model.normal_label: 1.0}


def gamete_distribution(
    model: InheritanceModel,
    sex: Sex,
    genotype: str,
) -> dict[str, float]:
    genotype = parse_genotype(model, sex, genotype)
    if model.mode == InheritanceMode.NON_GENETIC:
        return {"NA": 1.0}

    tokens = genotype.split("/")
    if model.mode in {
        InheritanceMode.Z_LINKED_RECESSIVE,
        InheritanceMode.Z_LINKED_DOMINANT,
    } and sex == Sex.FEMALE:
        return {tokens[0]: 0.5, "W": 0.5}

    counts = Counter(tokens)
    total = len(tokens)
    return {token: count / total for token, count in counts.items()}


def compatible_genotypes(
    model: InheritanceModel,
    sex: Sex,
    observed_phenotype: str,
) -> list[str]:
    """Enumerate parent genotypes compatible with a visible phenotype label.

    Carrier status is not inferred from appearance; a recessive heterozygote is
    therefore compatible with the normal visible label.
    """

    if observed_phenotype in {"", "unknown", "UNK", None}:  # type: ignore[arg-type]
        return enumerate_genotypes(model, sex)

    compatible: list[str] = []
    for genotype in enumerate_genotypes(model, sex):
        labels = phenotype_distribution(model, sex, genotype)
        if labels.get(observed_phenotype, 0.0) > 0.0:
            compatible.append(genotype)
    return compatible


def normalize_distribution(
    values: dict[str, float],
    *,
    allowed_keys: Iterable[str] | None = None,
) -> dict[str, float]:
    if not values:
        raise ValidationError("probability_distribution_empty")
    allowed = set(allowed_keys) if allowed_keys is not None else None
    normalized: dict[str, float] = {}
    total = 0.0
    for key, raw in values.items():
        if allowed is not None and key not in allowed:
            raise ValidationError(f"probability_distribution_key_invalid:{key}")
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ValidationError("probability_distribution_invalid")
        probability = float(raw)
        if not math.isfinite(probability):
            raise ValidationError("probability_distribution_non_finite")
        if probability < 0.0:
            raise ValidationError("probability_distribution_negative")
        total += probability
        normalized[key] = probability
    if total <= 0.0:
        raise ValidationError("probability_distribution_zero_total")
    return {key: probability / total for key, probability in normalized.items()}
