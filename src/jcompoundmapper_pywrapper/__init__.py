"""Wrapper for jCompoundMapper molecular fingerprints"""

from .jcompoundmapper import (
    DEFAULT_FP_PARAMETERS,
    AtomType,
    Fingerprint,
    Fingerprints_2D,
    Fingerprints_3D,
    FpParamType,
    JCompoundMapper,
)

__all__ = [
    "DEFAULT_FP_PARAMETERS",
    "AtomType",
    "Fingerprint",
    "Fingerprints_2D",
    "Fingerprints_3D",
    "FpParamType",
    "JCompoundMapper",
]

__version__ = "0.1.0"
