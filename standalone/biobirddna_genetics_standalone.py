#!/usr/bin/env python3
"""Single-file BioBirdDNA genetics prototype.

Use this file when a coding model or reviewer wants one portable script before
adopting the modular skill package. It supports:

- autosomal recessive/dominant/incomplete-dominant crosses;
- avian Z-linked recessive/dominant crosses under ZZ male / ZW female;
- uncertain parent genotype distributions;
- basic likelihood ranking of candidate inheritance models;
- a phenotype-only pigment rule screen that never confirms a gene.

Python standard library only.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import Enum
from itertools import product
from pathlib import Path
from typing import Any, Mapping


class Mode(str, Enum):
    AR = "autosomal_recessive"
    AD = "autosomal_dominant"
    INC = "incomplete_dominant"
    ZR = "z_linked_recessive"
    ZD = "z_linked_dominant"
    NON = "non_genetic"


@dataclass(frozen=True)
class Model:
    model_id: str
    mode: Mode
    normal: str = "N"
    variant: str = "m"
    normal_label: str = "normal"
    target_label: str = "target_colour"
    intermediate_label: str = "intermediate"
    penetrance: float = 1.0
    error: float = 0.01
    female_probability: float = 0.5
    baseline_target_probability: float = 0.05

    def __post_init__(self) -> None:
        if not 0.0 <= self.female_probability <= 1.0:
            raise ValueError("female_probability_must_be_between_0_and_1")
        # Z-linked offspring sex arises from maternal Z/W segregation, so the
        # simulated sex ratio is fixed at 0.5 and a configurable
        # female_probability cannot be honored. Fail closed instead of
        # silently ignoring the value (mirrors
        # src/biobirddna_genetics/models.py InheritanceModel.__post_init__
        # and npm/src/index.mjs modelFrom — keep in sync).
        if self.mode in {Mode.ZR, Mode.ZD} and self.female_probability != 0.5:
            raise ValueError("z_linked_female_probability_must_be_0_5")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "Model":
        return cls(
            model_id=str(value.get("model_id", value["mode"])),
            mode=Mode(value["mode"]),
            normal=str(value.get("normal_allele", "N")),
            variant=str(value.get("variant_allele", "m")),
            normal_label=str(value.get("normal_label", "normal")),
            target_label=str(value.get("target_label", "target_colour")),
            intermediate_label=str(value.get("intermediate_label", "intermediate")),
            penetrance=float(value.get("penetrance", 1.0)),
            error=float(value.get("observation_error", 0.01)),
            female_probability=float(value.get("female_probability", 0.5)),
            baseline_target_probability=float(
                value.get("baseline_target_probability", 0.05)
            ),
        )


def canonical(model: Model, genotype: str, sex: str) -> str:
    if model.mode == Mode.NON:
        return "NA"
    tokens = genotype.split("/")
    allowed = {model.normal, model.variant}
    if model.mode in {Mode.ZR, Mode.ZD} and sex == "female":
        if len(tokens) != 2 or tokens[1] != "W" or tokens[0] not in allowed:
            raise ValueError("z_linked_female_genotype_invalid")
        return genotype
    if len(tokens) != 2 or any(token not in allowed for token in tokens):
        raise ValueError("genotype_invalid")
    order = {model.normal: 0, model.variant: 1}
    return "/".join(sorted(tokens, key=lambda token: order[token]))


def gametes(model: Model, genotype: str, sex: str) -> dict[str, float]:
    genotype = canonical(model, genotype, sex)
    if genotype == "NA":
        return {"NA": 1.0}
    tokens = genotype.split("/")
    if model.mode in {Mode.ZR, Mode.ZD} and sex == "female":
        return {tokens[0]: 0.5, "W": 0.5}
    counts = Counter(tokens)
    return {token: count / 2.0 for token, count in counts.items()}


def phenotype(model: Model, genotype: str, sex: str) -> dict[str, float]:
    if model.mode == Mode.NON:
        p = model.baseline_target_probability
        return {model.target_label: p, model.normal_label: 1 - p}
    genotype = canonical(model, genotype, sex)
    copies = Counter(token for token in genotype.split("/") if token != "W")[
        model.variant
    ]
    if model.mode in {Mode.AR, Mode.ZR}:
        required = 1 if model.mode == Mode.ZR and sex == "female" else 2
        affected = copies >= required
    elif model.mode in {Mode.AD, Mode.ZD}:
        affected = copies >= 1
    elif model.mode == Mode.INC:
        if copies == 1:
            return {model.intermediate_label: 1.0}
        affected = copies == 2
    else:
        affected = False
    if affected:
        return {
            model.target_label: model.penetrance,
            model.normal_label: 1.0 - model.penetrance,
        }
    return {model.normal_label: 1.0}


def cross_known(model: Model, male: str, female: str) -> dict[str, Any]:
    if model.mode == Mode.NON:
        p = model.baseline_target_probability
        f = model.female_probability
        return {
            "phenotype_probabilities": {
                f"male:{model.target_label}": (1 - f) * p,
                f"male:{model.normal_label}": (1 - f) * (1 - p),
                f"female:{model.target_label}": f * p,
                f"female:{model.normal_label}": f * (1 - p),
            },
            "genotype_probabilities": {},
        }
    male = canonical(model, male, "male")
    female = canonical(model, female, "female")
    fg = gametes(model, male, "male")
    mg = gametes(model, female, "female")
    genotypes: dict[str, float] = defaultdict(float)
    phenotypes: dict[str, float] = defaultdict(float)

    for paternal, pp in fg.items():
        for maternal, mp in mg.items():
            base = pp * mp
            if model.mode in {Mode.ZR, Mode.ZD}:
                if maternal == "W":
                    sex = "female"
                    child = canonical(model, f"{paternal}/W", sex)
                else:
                    sex = "male"
                    child = canonical(model, f"{paternal}/{maternal}", sex)
                genotypes[f"{sex}:{child}"] += base
                for label, conditional in phenotype(model, child, sex).items():
                    phenotypes[f"{sex}:{label}"] += base * conditional
            else:
                child = canonical(model, f"{paternal}/{maternal}", "male")
                for sex, sex_probability in (
                    ("male", 1 - model.female_probability),
                    ("female", model.female_probability),
                ):
                    probability = base * sex_probability
                    genotypes[f"{sex}:{child}"] += probability
                    for label, conditional in phenotype(model, child, sex).items():
                        phenotypes[f"{sex}:{label}"] += probability * conditional
    return {
        "genotype_probabilities": dict(sorted(genotypes.items())),
        "phenotype_probabilities": dict(sorted(phenotypes.items())),
    }


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    total = sum(float(value) for value in values.values())
    if total <= 0:
        raise ValueError("distribution_zero_total")
    return {str(key): float(value) / total for key, value in values.items()}


def cross_uncertain(
    model: Model,
    males: Mapping[str, float],
    females: Mapping[str, float],
) -> dict[str, Any]:
    male_dist = normalize(males)
    female_dist = normalize(females)
    genotypes: dict[str, float] = defaultdict(float)
    phenotypes: dict[str, float] = defaultdict(float)
    for male, female in product(male_dist, female_dist):
        weight = male_dist[male] * female_dist[female]
        result = cross_known(model, male, female)
        for key, probability in result["genotype_probabilities"].items():
            genotypes[key] += weight * probability
        for key, probability in result["phenotype_probabilities"].items():
            phenotypes[key] += weight * probability
    return {
        "genotype_probabilities": dict(sorted(genotypes.items())),
        "phenotype_probabilities": dict(sorted(phenotypes.items())),
    }


def parent_genotypes(model: Model, sex: str, visible: str) -> list[str]:
    if model.mode == Mode.NON:
        return ["NA"]
    all_values = (
        [f"{model.normal}/W", f"{model.variant}/W"]
        if model.mode in {Mode.ZR, Mode.ZD} and sex == "female"
        else [
            f"{model.normal}/{model.normal}",
            f"{model.normal}/{model.variant}",
            f"{model.variant}/{model.variant}",
        ]
    )
    if visible in {"unknown", "UNK", ""}:
        return all_values
    return [
        value
        for value in all_values
        if phenotype(model, value, sex).get(visible, 0.0) > 0.0
    ]


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    if math.isinf(maximum):
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def infer(request: Mapping[str, Any]) -> dict[str, Any]:
    target = str(request.get("target_label", "target_colour"))
    models = [
        Model("M1_AR", Mode.AR, target_label=target),
        Model("M2_AD", Mode.AD, target_label=target),
        Model("M3_INC", Mode.INC, target_label=target),
        Model("M4_ZR", Mode.ZR, target_label=target),
        Model("M5_ZD", Mode.ZD, target_label=target),
        Model("M0_NONGENETIC", Mode.NON, target_label=target),
    ]
    counts = {
        f"{item.get('sex', 'unknown')}:{item['phenotype']}": int(item["count"])
        for item in request["offspring"]
    }
    sex_specific = all(not key.startswith("unknown:") for key in counts)
    rankings = []
    log_scores: dict[str, float] = {}

    for model in models:
        configs = []
        if model.mode == Mode.NON:
            configs = [("NA", "NA")]
        else:
            males = parent_genotypes(
                model,
                "male",
                str(request.get("male_parent_phenotype", "unknown")),
            )
            females = parent_genotypes(
                model,
                "female",
                str(request.get("female_parent_phenotype", "unknown")),
            )
            configs = list(product(males, females))
        terms = []
        config_rows = []
        for male, female in configs:
            result = cross_known(model, male, female)
            raw = result["phenotype_probabilities"]
            if not sex_specific:
                collapsed: dict[str, float] = defaultdict(float)
                for key, probability in raw.items():
                    collapsed[key.split(":", 1)[1]] += probability
                raw = {f"unknown:{key}": value for key, value in collapsed.items()}
            categories = sorted(set(counts) | set(raw))
            smoothed = {
                category: (1 - model.error) * raw.get(category, 0.0)
                + model.error / len(categories)
                for category in categories
            }
            total = sum(smoothed.values())
            smoothed = {key: value / total for key, value in smoothed.items()}
            log_likelihood = 0.0
            for category, count in counts.items():
                if count:
                    log_likelihood += count * math.log(smoothed[category])
            prior = 1.0 / max(1, len(configs))
            log_joint = math.log(prior) + log_likelihood
            terms.append(log_joint)
            config_rows.append(
                {
                    "male_genotype": male,
                    "female_genotype": female,
                    "log_likelihood": log_likelihood,
                    "predicted": smoothed,
                }
            )
        evidence = logsumexp(terms) if terms else float("-inf")
        log_scores[model.model_id] = evidence
        rankings.append(
            {
                "model_id": model.model_id,
                "mode": model.mode.value,
                "log_evidence": evidence,
                "configurations": config_rows,
            }
        )

    denominator = logsumexp(list(log_scores.values()))
    posterior = {
        model_id: math.exp(score - denominator)
        for model_id, score in log_scores.items()
    }
    for row in rankings:
        row["posterior"] = posterior[row["model_id"]]
    rankings.sort(key=lambda row: row["posterior"], reverse=True)
    return {
        "ranking": rankings,
        "best_model": rankings[0]["model_id"],
        "limits": {
            "claim": "candidate_inheritance_model_only",
            "gene_status": "not_assessed",
        },
    }


def pigment(request: Mapping[str, Any]) -> dict[str, Any]:
    limits = {
        "measurement_domain": "visible_phenotype_only",
        "biochemical_pigment_status": "unconfirmed",
        "gene_status": "not_assessed",
        "inheritance_status": "not_assessed",
    }
    image_quality = str(request.get("image_quality", "Q0"))
    if image_quality in {"Q0", "Q1"}:
        return {
            "primary_hypothesis": "UNK",
            "scores": {"UNK": 1.0},
            "rationale": {"UNK": ["image_quality_too_low_for_measurement"]},
            "needs_recapture": True,
            "limits": limits,
        }

    scores = defaultdict(lambda: 0.05)
    white = float(request.get("white_area_ratio", 0.0))
    if request.get("progressive_across_molts"):
        scores["MEL-PGR"] += 2.0
    if 0.02 <= white < 0.95:
        scores["MEL-LEU"] += 0.8 + white
    if white >= 0.95 and request.get("iris_red_confirmed") and not request.get(
        "flash_red_eye_risk"
    ):
        scores["MEL-ALB"] += 2.4
    if request.get("flash_red_eye_risk"):
        scores["IMG-ART"] += 1.4
    if request.get("localized_after_injury"):
        scores["LOC-ACQ"] += 1.8
    scores["MEL-BRN"] += 2.0 * float(request.get("brown_shift_score", 0.0))
    scores["MEL-DIL"] += 2.0 * float(request.get("grey_shift_score", 0.0))
    scores["MEL-MEL"] += 2.0 * float(request.get("dark_expansion_score", 0.0))
    scores["CAR-ALT"] += 2.0 * float(request.get("red_patch_change_score", 0.0))
    total = sum(scores.values())
    ranked = sorted(
        ((key, value / total) for key, value in scores.items()),
        key=lambda item: item[1],
        reverse=True,
    )
    return {
        "primary_hypothesis": ranked[0][0],
        "scores": dict(ranked),
        "needs_recapture": False,
        "limits": limits,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("operation", choices=["cross", "infer", "pigment"])
    parser.add_argument("--input", required=True, type=Path)
    args = parser.parse_args(argv)
    request = json.loads(args.input.read_text(encoding="utf-8"))
    if args.operation == "cross":
        model = Model.from_dict(request["model"])
        if "male_genotypes" in request:
            result = cross_uncertain(
                model,
                request["male_genotypes"],
                request["female_genotypes"],
            )
        else:
            result = cross_known(
                model,
                request["male_genotype"],
                request["female_genotype"],
            )
    elif args.operation == "infer":
        result = infer(request)
    else:
        result = pigment(request)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
