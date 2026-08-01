class BioBirdDNAError(Exception):
    """Base error for the BioBirdDNA genetics engine."""


class ValidationError(BioBirdDNAError):
    """Raised when a request or domain object violates an invariant."""


class UnsupportedModelError(BioBirdDNAError):
    """Raised when a model type is outside the implemented prototype."""


class EvidenceBoundaryError(BioBirdDNAError):
    """Raised when an output attempts a claim above its evidence level."""
