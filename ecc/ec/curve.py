# ecc/ec/curve.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..fields.base import Field, Element
from .point import Point


@dataclass(frozen=True)
class Curve:
    """
    Short Weierstrass curve: y^2 = x^3 + a x + b over an arbitrary field.
      - Over fields of characteristic 2 or 3, short Weierstrass form needs special care.
      - p > 3 is not enforced.
    """

    field: Field
    a: Element
    b: Element

    def __post_init__(self):
        F = self.field
        object.__setattr__(self, "a", F.coerce(self.a))
        object.__setattr__(self, "b", F.coerce(self.b))

        # Nonsingularity check when meaningful:
        # Note that short Weierstrass assumes char(F) != 2 or 3
        p = getattr(F, "p", None)
        if isinstance(p, int) and p in (2, 3):
            raise ValueError(
                "Short Weierstrass form y^2 = x^3 + ax + b with the standard addition formulas "
                "is not supported in characteristic 2 or 3. "
                "Use a different curve model for p=2 or p=3."
            )
        # discriminant ∝ 4a^3 + 27b^2 must be nonzero
        disc = (F.coerce(4) * (self.a ** 3)) + (F.coerce(27) * (self.b ** 2))
        if disc == F.zero:
            raise ValueError("Singular curve (discriminant is zero)")

    def infinity(self) -> Point:
        return Point.infinity(self)

    def point(self, x: Any, y: Any, *, check: bool = True) -> Point:
        P = Point(self, x, y, is_infinity=False)
        if check and not self.is_on_curve(P):
            raise ValueError("Point is not on this curve")
        return P

    def is_on_curve(self, P: Point) -> bool:
        if P.is_infinity:
            return True
        x, y = P.x, P.y
        assert x is not None and y is not None
        left = y ** 2
        right = (x ** 3) + (self.a * x) + self.b
        return left == right

    def neg(self, P: Point) -> Point:
        if P.is_infinity:
            return P
        assert P.x is not None and P.y is not None
        return Point(self, P.x, -P.y)

    def add(self, P: Point, Q: Point) -> Point:
        if P.curve is not self or Q.curve is not self:
            raise TypeError("Cannot add points from different curves")

        F = self.field

        if P.is_infinity:
            return Q
        if Q.is_infinity:
            return P

        x1, y1 = P.x, P.y
        x2, y2 = Q.x, Q.y
        assert x1 is not None and y1 is not None and x2 is not None and y2 is not None

        # P + (-P) = ∞  (vertical line)
        if x1 == x2 and (y1 + y2) == F.zero:
            return self.infinity()

        if x1 != x2:
            lam = (y2 - y1) / (x2 - x1)
        else:
            # Doubling: if y == 0, tangent is vertical -> infinity
            if y1 == F.zero:
                return self.infinity()
            lam = (F.coerce(3) * (x1 ** 2) + self.a) / (F.coerce(2) * y1)

        x3 = (lam ** 2) - x1 - x2
        y3 = lam * (x1 - x3) - y1
        return Point(self, x3, y3)

    def scalar_mul(self, k: int, P: Point) -> Point:
        if P.curve is not self:
            raise TypeError("Point is not on this curve")
        if k == 0 or P.is_infinity:
            return self.infinity()
        if k < 0:
            return self.scalar_mul(-k, self.neg(P))

        result = self.infinity()
        addend = P
        n = k
        while n:
            if n & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
            n >>= 1
        return result
    