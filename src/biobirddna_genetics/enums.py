from __future__ import annotations

from enum import Enum


class Sex(str, Enum):
    MALE = "male"
    FEMALE = "female"
    UNKNOWN = "unknown"


class InheritanceMode(str, Enum):
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    INCOMPLETE_DOMINANT = "incomplete_dominant"
    Z_LINKED_RECESSIVE = "z_linked_recessive"
    Z_LINKED_DOMINANT = "z_linked_dominant"
    NON_GENETIC = "non_genetic"


class GenotypeState(str, Enum):
    NORMAL = "normal"
    CARRIER = "carrier"
    INTERMEDIATE = "intermediate"
    AFFECTED = "affected"
    UNKNOWN = "unknown"


class EvidenceGrade(str, Enum):
    G0 = "G0"  # image/appearance only
    G1 = "G1"  # one pair or one brood
    G2 = "G2"  # repeated broods from the same pair
    G3 = "G3"  # independent pairs and reciprocal crosses
    G4 = "G4"  # multi-generation pedigree and test cross
    G5 = "G5"  # molecular cosegregation
    G6 = "G6"  # independent replication or functional evidence


class PhenotypeEvidence(str, Enum):
    E0 = "E0"
    E1 = "E1"
    E2 = "E2"
    E3 = "E3"
    E4 = "E4"
    E5 = "E5"


class ImageQuality(str, Enum):
    Q0 = "Q0"
    Q1 = "Q1"
    Q2 = "Q2"
    Q3 = "Q3"
    Q4 = "Q4"


class ConfidenceLevel(str, Enum):
    C0 = "C0"
    C1 = "C1"
    C2 = "C2"
    C3 = "C3"
    C4 = "C4"
