# ecc/ec/point.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from ..fields.base import Element


@dataclass(frozen=True)
class Point:
    """
    Affine point on an elliptic curve, or the point at infinity.
      - if is_infinity: x=y=None
      - else: x,y are Elements in curve.field
    """

    curve: Any
    x: Optional[Element] = None
    y: Optional[Element] = None
    is_infinity: bool = False

    @staticmethod
    def infinity(curve: Any) -> "Point":
        return Point(curve=curve, x=None, y=None, is_infinity=True)

    def __post_init__(self):
        if self.is_infinity:
            if self.x is not None or self.y is not None:
                raise ValueError("Infinity point must have x=y=None")
            return

        if self.x is None or self.y is None:
            raise ValueError("Non-infinity point must have x and y")

        F = self.curve.field
        # Coerce into the curve's field
        object.__setattr__(self, "x", F.coerce(self.x))
        object.__setattr__(self, "y", F.coerce(self.y))

    def __neg__(self) -> "Point":
        return self.curve.neg(self)

    def __add__(self, other: "Point") -> "Point":
        return self.curve.add(self, other)

    def __radd__(self, other: "Point") -> "Point":
        return self.__add__(other)

    def __sub__(self, other: "Point") -> "Point":
        return self + (-other)

    def __mul__(self, k: int) -> "Point":
        if not isinstance(k, int):
            raise TypeError("Can only multiply a Point by an int scalar")
        return self.curve.scalar_mul(k, self)

    def __rmul__(self, k: int) -> "Point":
        return self.__mul__(k)

    def __repr__(self) -> str:
        if self.is_infinity:
            return "Point(∞)"
        return f"Point({self.x}, {self.y})"