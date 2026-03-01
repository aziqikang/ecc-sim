# ecc/__init__.py
"""
ecc: minimal elliptic curve playground over interchangeable fields

Public API:
- Fields: RationalField, PrimeField
- EC: Curve, Point
"""

from .fields import RationalField, PrimeField
from .ec import Curve, Point

__all__ = ["RationalField", "PrimeField", "Curve", "Point"]