"""BioBirdDNA Genetics deterministic and probabilistic prototype.

The package intentionally separates visible phenotype, inheritance-model inference,
and molecular confirmation. It must not promote an image-only observation into a
confirmed gene claim.
"""

from .models import InheritanceModel
from .breeding_records import normalize_breeding_records
from .simulation import simulate_recovery
from .cross import cross_known, cross_uncertain
from .inference import infer_candidate_models
from .pedigree import Pedigree
from .optimizer import rank_mating_pairs
from .pigment import rank_pigment_hypotheses

__all__ = [
    "InheritanceModel",
    "normalize_breeding_records",
    "simulate_recovery",
    "cross_known",
    "cross_uncertain",
    "infer_candidate_models",
    "Pedigree",
    "rank_mating_pairs",
    "rank_pigment_hypotheses",
]

__version__ = "0.2.0"
