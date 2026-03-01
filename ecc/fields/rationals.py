# ecc/fields/rationals.py
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .base import Element


class RationalField:
    """The field of rational numbers Q, backed by fractions.Fraction."""

    def __init__(self):
        # singletons for convenience
        self._zero = Rational(self, Fraction(0))
        self._one = Rational(self, Fraction(1))

    @property
    def zero(self) -> "Rational":
        return self._zero

    @property
    def one(self) -> "Rational":
        return self._one

    def coerce(self, x: Any) -> "Rational":
        if isinstance(x, Rational):
            if x.field is not self:
                raise TypeError("Field mismatch: Rational from a different RationalField context")
            return x
        if isinstance(x, Fraction):
            return Rational(self, x)
        if isinstance(x, int):
            return Rational(self, Fraction(x))
        # allow "int-like" strings? keep strict for now
        raise TypeError(f"Cannot coerce type {type(x).__name__} into RationalField")

    def __call__(self, x: Any) -> "Rational":
        return self.coerce(x)

    def __repr__(self) -> str:
        return "RationalField(Q)"


@dataclass(frozen=True)
class Rational(Element):
    field: RationalField
    value: Fraction

    def _add(self, other: "Rational") -> "Rational":
        return Rational(self.field, self.value + other.value)

    def _sub(self, other: "Rational") -> "Rational":
        return Rational(self.field, self.value - other.value)

    def _mul(self, other: "Rational") -> "Rational":
        return Rational(self.field, self.value * other.value)

    def _neg(self) -> "Rational":
        return Rational(self.field, -self.value)

    def _inv(self) -> "Rational":
        if self.value == 0:
            raise ZeroDivisionError("inverse of 0 does not exist in Q")
        return Rational(self.field, 1 / self.value)

    def __repr__(self) -> str:
        return f"Q({self.value})"
    