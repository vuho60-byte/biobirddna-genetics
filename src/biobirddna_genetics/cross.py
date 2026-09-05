from __future__ import annotations

from collections import defaultdict
from typing import Any, Mapping

from .enums import InheritanceMode, Sex
from .errors import ValidationError
from .genotype import (
    gamete_distribution,
    normalize_distribution,
    parse_genotype,
    phenotype_distribution,
)
from .models import InheritanceModel


def _add_phenotype(
    destination: dict[str, float],
    model: InheritanceModel,
    sex: Sex,
    genotype: str,
    probability: float,
) -> None:
    for phenotype, conditional_probability in phenotype_distribution(
        model,
        sex,
        genotype,
    ).items():
        destination[f"{sex.value}:{phenotype}"] += (
            probability * conditional_probability
        )


def cross_known(
    model: InheritanceModel,
    male_genotype: str,
    female_genotype: str,
) -> dict[str, Any]:
    """Simulate one cross when both parent genotypes are known."""

    male_genotype = parse_genotype(model, Sex.MALE, male_genotype)
    female_genotype = parse_genotype(model, Sex.FEMALE, female_genotype)

    if model.mode == InheritanceMode.NON_GENETIC:
        p_female = model.female_probability
        phenotype_probabilities = {
            f"male:{model.target_label}": (1.0 - p_female)
            * model.baseline_target_probability,
            f"male:{model.normal_label}": (1.0 - p_female)
            * (1.0 - model.baseline_target_probability),
            f"female:{model.target_label}": p_female
            * model.baseline_target_probability,
            f"female:{model.normal_label}": p_female
            * (1.0 - model.baseline_target_probability),
        }
        return {
            "model_id": model.model_id,
            "mode": model.mode.value,
            "male_genotype": "NA",
            "female_genotype": "NA",
            "genotype_probabilities": {},
            "phenotype_probabilities": phenotype_probabilities,
            "limits": {
                "gene_status": "not_assessed",
                "interpretation": "non_genetic_baseline_model",
            },
        }

    father_gametes = gamete_distribution(model, Sex.MALE, male_genotype)
    mother_gametes = gamete_distribution(model, Sex.FEMALE, female_genotype)
    genotype_probabilities: dict[str, float] = defaultdict(float)
    phenotype_probabilities: dict[str, float] = defaultdict(float)

    if model.mode in {
        InheritanceMode.AUTOSOMAL_RECESSIVE,
        InheritanceMode.AUTOSOMAL_DOMINANT,
        InheritanceMode.INCOMPLETE_DOMINANT,
    }:
        p_female = model.female_probability
        for paternal_allele, paternal_probability in father_gametes.items():
            for maternal_allele, maternal_probability in mother_gametes.items():
                base_probability = paternal_probability * maternal_probability
                genotype = parse_genotype(
                    model,
                    Sex.MALE,
                    f"{paternal_allele}/{maternal_allele}",
                )
                for sex, sex_probability in (
                    (Sex.MALE, 1.0 - p_female),
                    (Sex.FEMALE, p_female),
                ):
                    probability = base_probability * sex_probability
                    genotype_probabilities[f"{sex.value}:{genotype}"] += probability
                    _add_phenotype(
                        phenotype_probabilities,
                        model,
                        sex,
                        genotype,
                        probability,
                    )

    elif model.mode in {
        InheritanceMode.Z_LINKED_RECESSIVE,
        InheritanceMode.Z_LINKED_DOMINANT,
    }:
        for paternal_z, paternal_probability in father_gametes.items():
            for maternal_chromosome, maternal_probability in mother_gametes.items():
                probability = paternal_probability * maternal_probability
                if maternal_chromosome == "W":
                    sex = Sex.FEMALE
                    genotype = parse_genotype(model, sex, f"{paternal_z}/W")
                else:
                    sex = Sex.MALE
                    genotype = parse_genotype(
                        model,
                        sex,
                        f"{paternal_z}/{maternal_chromosome}",
                    )
                genotype_probabilities[f"{sex.value}:{genotype}"] += probability
                _add_phenotype(
                    phenotype_probabilities,
                    model,
                    sex,
                    genotype,
                    probability,
                )
    else:
        raise ValidationError(f"cross_mode_not_supported:{model.mode.value}")

    return {
        "model_id": model.model_id,
        "mode": model.mode.value,
        "male_genotype": male_genotype,
        "female_genotype": female_genotype,
        "genotype_probabilities": dict(sorted(genotype_probabilities.items())),
        "phenotype_probabilities": dict(sorted(phenotype_probabilities.items())),
        "limits": {
            "gene_status": "model_assumption_only",
            "inheritance_status": "simulation_not_confirmation",
        },
    }


def cross_uncertain(
    model: InheritanceModel,
    male_genotypes: Mapping[str, float],
    female_genotypes: Mapping[str, float],
) -> dict[str, Any]:
    """Marginalize a cross over uncertain parent genotype distributions."""

    parsed_male: dict[str, float] = defaultdict(float)
    parsed_female: dict[str, float] = defaultdict(float)
    for genotype, probability in normalize_distribution(dict(male_genotypes)).items():
        parsed_male[parse_genotype(model, Sex.MALE, genotype)] += probability
    for genotype, probability in normalize_distribution(dict(female_genotypes)).items():
        parsed_female[parse_genotype(model, Sex.FEMALE, genotype)] += probability
    male_distribution = normalize_distribution(dict(parsed_male))
    female_distribution = normalize_distribution(dict(parsed_female))

    if model.mode == InheritanceMode.NON_GENETIC:
        return cross_known(model, "NA", "NA")

    genotype_probabilities: dict[str, float] = defaultdict(float)
    phenotype_probabilities: dict[str, float] = defaultdict(float)
    component_crosses: list[dict[str, Any]] = []

    for male_genotype, male_probability in male_distribution.items():
        for female_genotype, female_probability in female_distribution.items():
            pair_probability = male_probability * female_probability
            result = cross_known(model, male_genotype, female_genotype)
            component_crosses.append(
                {
                    "male_genotype": male_genotype,
                    "female_genotype": female_genotype,
                    "weight": pair_probability,
                }
            )
            for key, probability in result["genotype_probabilities"].items():
                genotype_probabilities[key] += pair_probability * probability
            for key, probability in result["phenotype_probabilities"].items():
                phenotype_probabilities[key] += pair_probability * probability

    return {
        "model_id": model.model_id,
        "mode": model.mode.value,
        "male_genotype_posterior": male_distribution,
        "female_genotype_posterior": female_distribution,
        "genotype_probabilities": dict(sorted(genotype_probabilities.items())),
        "phenotype_probabilities": dict(sorted(phenotype_probabilities.items())),
        "component_crosses": component_crosses,
        "limits": {
            "gene_status": "inferred_parent_genotypes",
            "inheritance_status": "simulation_not_confirmation",
        },
    }


def collapse_phenotypes_by_sex(
    phenotype_probabilities: Mapping[str, float],
) -> dict[str, float]:
    collapsed: dict[str, float] = defaultdict(float)
    for key, probability in phenotype_probabilities.items():
        try:
            _sex, phenotype = key.split(":", 1)
        except ValueError as exc:
            raise ValidationError(f"phenotype_probability_key_invalid:{key}") from exc
        collapsed[phenotype] += float(probability)
    return dict(sorted(collapsed.items()))
