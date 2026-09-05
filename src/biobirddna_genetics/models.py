from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any, Mapping

from .enums import InheritanceMode, Sex
from .errors import ValidationError


def _probability(value: float, field_name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValidationError(f"{field_name}_must_be_between_0_and_1")
    return value


def _nonnegative_finite(value: float, field_name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or float(value) < 0.0
    ):
        raise ValidationError(f"{field_name}_invalid")
    return float(value)


@dataclass(frozen=True)
class InheritanceModel:
    model_id: str
    mode: InheritanceMode
    normal_allele: str = "N"
    variant_allele: str = "m"
    normal_label: str = "normal"
    target_label: str = "target_colour"
    intermediate_label: str = "intermediate"
    penetrance: float = 1.0
    observation_error: float = 0.01
    female_probability: float = 0.5
    baseline_target_probability: float = 0.05
    description: str = ""

    def __post_init__(self) -> None:
        if not self.model_id.strip():
            raise ValidationError("model_id_required")
        if self.normal_allele == self.variant_allele:
            raise ValidationError("alleles_must_differ")
        phenotype_labels = {self.normal_label, self.target_label}
        if len(phenotype_labels) != 2 or (
            self.mode == InheritanceMode.INCOMPLETE_DOMINANT
            and self.intermediate_label in phenotype_labels
        ):
            raise ValidationError("phenotype_labels_must_be_distinct")
        object.__setattr__(self, "penetrance", _probability(self.penetrance, "penetrance"))
        object.__setattr__(
            self,
            "observation_error",
            _probability(self.observation_error, "observation_error"),
        )
        object.__setattr__(
            self,
            "female_probability",
            _probability(self.female_probability, "female_probability"),
        )
        object.__setattr__(
            self,
            "baseline_target_probability",
            _probability(
                self.baseline_target_probability,
                "baseline_target_probability",
            ),
        )
        # Z-linked offspring sex arises from maternal Z/W segregation, so the
        # simulated sex ratio is fixed at 0.5 and a configurable
        # female_probability cannot be honored. Fail closed instead of
        # silently ignoring the value (mirrored in npm/src/index.mjs
        # modelFrom — keep in sync).
        if (
            self.mode
            in {
                InheritanceMode.Z_LINKED_RECESSIVE,
                InheritanceMode.Z_LINKED_DOMINANT,
            }
            and self.female_probability != 0.5
        ):
            raise ValidationError("z_linked_female_probability_must_be_0_5")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "InheritanceModel":
        allowed = {
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
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise ValidationError(f"model_unknown_field:{unknown[0]}")

        def required_text(name: str, default: str | None = None) -> str:
            raw = value.get(name, default)
            if (
                not isinstance(raw, str)
                or not raw
                or raw != raw.strip()
            ):
                raise ValidationError(f"{name}_required")
            return raw

        def probability(name: str, default: float) -> float:
            raw = value.get(name, default)
            if (
                isinstance(raw, bool)
                or not isinstance(raw, (int, float))
                or not math.isfinite(float(raw))
                or not 0.0 <= float(raw) <= 1.0
            ):
                raise ValidationError(f"{name}_invalid")
            return float(raw)

        try:
            mode = InheritanceMode(value["mode"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError("mode_invalid") from exc
        description = value.get("description", "")
        if not isinstance(description, str):
            raise ValidationError("description_invalid")
        return cls(
            model_id=required_text("model_id"),
            mode=mode,
            normal_allele=required_text("normal_allele", "N"),
            variant_allele=required_text("variant_allele", "m"),
            normal_label=required_text("normal_label", "normal"),
            target_label=required_text("target_label", "target_colour"),
            intermediate_label=required_text(
                "intermediate_label",
                "intermediate",
            ),
            penetrance=probability("penetrance", 1.0),
            observation_error=probability("observation_error", 0.01),
            female_probability=probability("female_probability", 0.5),
            baseline_target_probability=probability(
                "baseline_target_probability",
                0.05,
            ),
            description=description,
        )

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["mode"] = self.mode.value
        return value


@dataclass(frozen=True)
class OffspringObservation:
    phenotype: str
    count: int
    sex: Sex = Sex.UNKNOWN

    def __post_init__(self) -> None:
        if not self.phenotype.strip():
            raise ValidationError("offspring_phenotype_required")
        if self.count < 0:
            raise ValidationError("offspring_count_must_be_nonnegative")


@dataclass(frozen=True)
class PedigreeIndividual:
    individual_id: str
    sex: Sex = Sex.UNKNOWN
    sire_id: str | None = None
    dam_id: str | None = None

    def __post_init__(self) -> None:
        if not self.individual_id.strip():
            raise ValidationError("individual_id_required")
        if self.sire_id == self.individual_id or self.dam_id == self.individual_id:
            raise ValidationError("individual_cannot_be_own_parent")


@dataclass
class GeneticProfile:
    individual_id: str
    sex: Sex
    genotype_posteriors: dict[str, dict[str, float]] = field(default_factory=dict)
    phenotype_label: str = "unknown"
    health_risk: float = 0.0
    fertility_risk: float = 0.0
    welfare_blocked: bool = False
    notes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.individual_id.strip():
            raise ValidationError("individual_id_required")
        self.health_risk = _probability(self.health_risk, "health_risk")
        self.fertility_risk = _probability(self.fertility_risk, "fertility_risk")


@dataclass(frozen=True)
class ModelHypothesis:
    model: InheritanceModel
    weight: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "weight",
            _nonnegative_finite(self.weight, "model_weight"),
        )


@dataclass(frozen=True)
class PigmentObservation:
    white_area_ratio: float = 0.0
    brown_shift_score: float = 0.0
    grey_shift_score: float = 0.0
    dark_expansion_score: float = 0.0
    red_patch_change_score: float = 0.0
    pattern_retained: bool | None = None
    iris_red_confirmed: bool = False
    flash_red_eye_risk: bool = False
    progressive_across_molts: bool = False
    present_in_first_plumage: bool | None = None
    localized_after_injury: bool = False
    bilateral_symmetry: float | None = None
    image_quality: str = "Q0"

    def __post_init__(self) -> None:
        for field_name in (
            "white_area_ratio",
            "brown_shift_score",
            "grey_shift_score",
            "dark_expansion_score",
            "red_patch_change_score",
        ):
            _probability(getattr(self, field_name), field_name)
        if self.bilateral_symmetry is not None:
            _probability(self.bilateral_symmetry, "bilateral_symmetry")
        for field_name in (
            "iris_red_confirmed",
            "flash_red_eye_risk",
            "progressive_across_molts",
            "localized_after_injury",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValidationError(f"{field_name}_invalid")
        for field_name in (
            "pattern_retained",
            "present_in_first_plumage",
        ):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, bool):
                raise ValidationError(f"{field_name}_invalid")
        if self.image_quality not in {"Q0", "Q1", "Q2", "Q3", "Q4"}:
            raise ValidationError("image_quality_invalid")


@dataclass(frozen=True)
class PairingWeights:
    target: float = 1.0
    information: float = 0.7
    inbreeding: float = 2.0
    health: float = 1.5
    uncertainty: float = 0.4

    def __post_init__(self) -> None:
        for field_name in (
            "target",
            "information",
            "inbreeding",
            "health",
            "uncertainty",
        ):
            object.__setattr__(
                self,
                field_name,
                _nonnegative_finite(
                    getattr(self, field_name),
                    f"{field_name}_weight",
                ),
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any] | None) -> "PairingWeights":
        if not value:
            return cls()
        return cls(
            target=value.get("target", 1.0),
            information=value.get("information", 0.7),
            inbreeding=value.get("inbreeding", 2.0),
            health=value.get("health", 1.5),
            uncertainty=value.get("uncertainty", 0.4),
        )
