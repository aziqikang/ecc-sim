# ecc/fields/prime.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .base import Element


def _is_prime_trial(n: int) -> bool:
    """Simple deterministic primality check suitable for small/medium n."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


class PrimeField:
    """Prime field F_p, elements are integers mod p."""

    def __init__(self, p: int, *, validate_prime: bool = True):
        if not isinstance(p, int):
            raise TypeError("p must be an int")
        if p < 2:
            raise ValueError("p must be >= 2")
        if validate_prime and not _is_prime_trial(p):
            raise ValueError(f"p must be prime; got {p}")
        self.p = p
        self._zero = Fp(self, 0)
        self._one = Fp(self, 1 % p)

    @property
    def zero(self) -> "Fp":
        return self._zero

    @property
    def one(self) -> "Fp":
        return self._one

    def normalize(self, x: int) -> int:
        return x % self.p

    def coerce(self, x: Any) -> "Fp":
        if isinstance(x, Fp):
            if x.field is not self:
                raise TypeError("Field mismatch: Fp element from a different PrimeField context")
            return x
        if isinstance(x, int):
            return Fp(self, x % self.p)
        raise TypeError(f"Cannot coerce type {type(x).__name__} into PrimeField(F_{self.p})")

    def __call__(self, x: Any) -> "Fp":
        return self.coerce(x)

    def __repr__(self) -> str:
        return f"PrimeField(F_{self.p})"


@dataclass(frozen=True)
class Fp(Element):
    field: PrimeField
    value: int  # always reduced mod p

    def _add(self, other: "Fp") -> "Fp":
        p = self.field.p
        return Fp(self.field, (self.value + other.value) % p)

    def _sub(self, other: "Fp") -> "Fp":
        p = self.field.p
        return Fp(self.field, (self.value - other.value) % p)

    def _mul(self, other: "Fp") -> "Fp":
        p = self.field.p
        return Fp(self.field, (self.value * other.value) % p)

    def _neg(self) -> "Fp":
        p = self.field.p
        return Fp(self.field, (-self.value) % p)

    def _inv(self) -> "Fp":
        a = self.value
        p = self.field.p
        if a == 0:
            raise ZeroDivisionError("inverse of 0 does not exist in F_p")
        # Extended Euclidean Algorithm for a^{-1} mod p
        t, new_t = 0, 1
        r, new_r = p, a
        while new_r != 0:
            q = r // new_r
            t, new_t = new_t, t - q * new_t
            r, new_r = new_r, r - q * new_r
        # If p is prime and a != 0, r should be 1
        if r != 1:
            raise ZeroDivisionError("no inverse exists (is p prime?)")
        return Fp(self.field, t % p)

    def __repr__(self) -> str:
        return f"Fp_{self.field.p}({self.value})"
    