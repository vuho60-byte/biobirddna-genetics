from __future__ import annotations

from collections import defaultdict
from typing import Any

from .models import PigmentObservation


LABELS = (
    "MEL-LEU",
    "MEL-PGR",
    "MEL-ALB",
    "MEL-INO",
    "MEL-BRN",
    "MEL-DIL",
    "MEL-MEL",
    "CAR-ALT",
    "LOC-ACQ",
    "IMG-ART",
    "BASE-NORM",
    "UNK",
)


def _normalize(scores: dict[str, float]) -> dict[str, float]:
    clipped = {label: max(0.0, float(scores.get(label, 0.0))) for label in LABELS}
    total = sum(clipped.values())
    if total <= 0.0:
        return {label: (1.0 if label == "UNK" else 0.0) for label in LABELS}
    return {label: score / total for label, score in clipped.items()}


def rank_pigment_hypotheses(observation: PigmentObservation) -> dict[str, Any]:
    """Rule-based phenotype ranking only.

    This function ranks T02 working labels from visible features. It deliberately
    does not identify biochemical pigment molecules, genes, dominance, or mode of
    inheritance.
    """

    limits = {
        "measurement_domain": "visible_phenotype_only",
        "biochemical_pigment_status": "unconfirmed",
        "gene_status": "not_assessed",
        "inheritance_status": "not_assessed",
        "required_language": "kieu_hinh_goi_y",
    }
    if observation.image_quality in {"Q0", "Q1"}:
        return {
            "primary_hypothesis": "UNK",
            "alternative_hypotheses": [],
            "scores": _normalize({"UNK": 1.0}),
            "rationale": {"UNK": ["image_quality_too_low_for_measurement"]},
            "needs_recapture": True,
            "limits": limits,
        }

    scores: dict[str, float] = defaultdict(lambda: 0.05)
    rationale: dict[str, list[str]] = defaultdict(list)

    if observation.flash_red_eye_risk:
        scores["IMG-ART"] += 1.4
        rationale["IMG-ART"].append("flash_red_eye_risk")

    if observation.localized_after_injury:
        scores["LOC-ACQ"] += 1.8
        rationale["LOC-ACQ"].append("localized_after_injury")

    if observation.progressive_across_molts:
        scores["MEL-PGR"] += 2.0
        rationale["MEL-PGR"].append("white_or_pale_area_progresses_across_molts")

    if observation.white_area_ratio >= 0.95:
        if observation.iris_red_confirmed and not observation.flash_red_eye_risk:
            scores["MEL-ALB"] += 2.4
            rationale["MEL-ALB"].append("near_complete_white_with_confirmed_red_iris")
        elif observation.pattern_retained:
            scores["MEL-INO"] += 1.5
            rationale["MEL-INO"].append("very_pale_with_residual_pattern")
        else:
            scores["MEL-LEU"] += 1.0
            rationale["MEL-LEU"].append("near_complete_white_without_red_iris_confirmation")

    if 0.02 <= observation.white_area_ratio < 0.95:
        scores["MEL-LEU"] += 0.8 + 1.2 * observation.white_area_ratio
        rationale["MEL-LEU"].append("partial_white_plumage")
        if observation.present_in_first_plumage is False:
            scores["MEL-PGR"] += 0.7
            rationale["MEL-PGR"].append("not_present_in_first_plumage")

    if observation.brown_shift_score > 0.15:
        scores["MEL-BRN"] += 0.5 + 2.0 * observation.brown_shift_score
        rationale["MEL-BRN"].append("black_to_brown_visible_shift")
        if observation.pattern_retained:
            scores["MEL-BRN"] += 0.4
            rationale["MEL-BRN"].append("pattern_retained")

    if observation.grey_shift_score > 0.15:
        scores["MEL-DIL"] += 0.5 + 2.0 * observation.grey_shift_score
        rationale["MEL-DIL"].append("black_to_grey_or_silver_visible_shift")
        if observation.pattern_retained:
            scores["MEL-DIL"] += 0.4
            rationale["MEL-DIL"].append("pattern_retained")

    if observation.dark_expansion_score > 0.15:
        scores["MEL-MEL"] += 0.5 + 2.0 * observation.dark_expansion_score
        rationale["MEL-MEL"].append("dark_area_expansion")

    if observation.red_patch_change_score > 0.15:
        scores["CAR-ALT"] += 0.5 + 2.0 * observation.red_patch_change_score
        rationale["CAR-ALT"].append("red_or_orange_domain_changed")

    visible_signal = max(
        observation.white_area_ratio,
        observation.brown_shift_score,
        observation.grey_shift_score,
        observation.dark_expansion_score,
        observation.red_patch_change_score,
    )
    if visible_signal < 0.05 and not observation.localized_after_injury:
        scores["BASE-NORM"] += 1.5
        rationale["BASE-NORM"].append("no_strong_abnormal_visible_signal")

    if observation.image_quality == "Q2":
        scores["UNK"] += 0.5
        rationale["UNK"].append("image_quality_supports_visual_screening_only")

    normalized = _normalize(scores)
    ranking = sorted(normalized.items(), key=lambda item: item[1], reverse=True)
    primary = ranking[0][0]
    # Only labels with score > 0.0 are real candidates: the scoring rules only
    # produce exactly 0.0 when no rule fired for a label. Strict float
    # comparison is intentional — genuinely tiny positive scores are kept.
    alternatives = [label for label, score in ranking[1:4] if score > 0.0]

    return {
        "primary_hypothesis": primary,
        "alternative_hypotheses": alternatives,
        "scores": dict(ranking),
        "rationale": {label: reasons for label, reasons in rationale.items()},
        "needs_recapture": False,
        "limits": limits,
    }
