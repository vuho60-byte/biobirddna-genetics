from __future__ import annotations

from .enums import EvidenceGrade
from .errors import EvidenceBoundaryError


_EVIDENCE_ORDER = {
    EvidenceGrade.G0: 0,
    EvidenceGrade.G1: 1,
    EvidenceGrade.G2: 2,
    EvidenceGrade.G3: 3,
    EvidenceGrade.G4: 4,
    EvidenceGrade.G5: 5,
    EvidenceGrade.G6: 6,
}


_MINIMUM = {
    # Observation- and phenotype-level claims (scientific-boundaries.md
    # section 1, rows "Quan sát" and "Kiểu hình"): recording submitted
    # observations and suggesting phenotype working labels make no genetic
    # claim, so their floor is G0 (image/appearance only).
    "observation_record": EvidenceGrade.G0,
    "phenotype_hypothesis": EvidenceGrade.G0,
    "simulate_cross": EvidenceGrade.G0,
    "candidate_inheritance_model": EvidenceGrade.G1,
    "supported_inheritance_pattern": EvidenceGrade.G3,
    "strong_inheritance_model": EvidenceGrade.G4,
    "molecular_locus_supported": EvidenceGrade.G5,
    "mechanism_strongly_confirmed": EvidenceGrade.G6,
}


def claim_allowed(claim_type: str, grade: EvidenceGrade) -> bool:
    required = _MINIMUM.get(claim_type)
    if required is None:
        return False
    return _EVIDENCE_ORDER[grade] >= _EVIDENCE_ORDER[required]


def assert_claim_allowed(claim_type: str, grade: EvidenceGrade) -> None:
    if not claim_allowed(claim_type, grade):
        required = _MINIMUM.get(claim_type)
        required_text = required.value if required else "unsupported"
        raise EvidenceBoundaryError(
            f"claim_not_allowed:{claim_type}:requires_{required_text}:has_{grade.value}"
        )


def safe_claim_language(claim_type: str, grade: EvidenceGrade) -> str:
    if claim_allowed(claim_type, grade):
        return claim_type
    if _EVIDENCE_ORDER[grade] <= _EVIDENCE_ORDER[EvidenceGrade.G2]:
        return "exploratory_candidate_model"
    if _EVIDENCE_ORDER[grade] == _EVIDENCE_ORDER[EvidenceGrade.G3]:
        return "segregation_pattern_supported"
    if _EVIDENCE_ORDER[grade] == _EVIDENCE_ORDER[EvidenceGrade.G4]:
        return "strong_pedigree_model_without_molecular_confirmation"
    return "insufficient_evidence"
