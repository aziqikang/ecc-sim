# ecc/fields/base.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class Field(Protocol):
    """A field context that can coerce Python values into Elements of this field."""

    @property
    def zero(self) -> "Element": ...
    @property
    def one(self) -> "Element": ...

    def coerce(self, x: Any) -> "Element": ...

    # Convenience: allow calling the field like F(x)
    def __call__(self, x: Any) -> "Element": ...


@dataclass(frozen=True)
class Element:
    """
    Generic field element base class.

    Concrete subclasses should:
      - store a 'field' reference
      - store a canonical 'value' representation
      - implement _add/_sub/_mul/_neg/_inv/_eq/_hash consistently
    """

    field: Any  # concrete Field type in subclasses
    value: Any

    # ---- coercion helpers ----
    def _coerce_other(self, other: Any) -> "Element":
        return self.field.coerce(other)

    def _check_same_field(self, other: "Element") -> None:
        if self.field is not other.field:
            raise TypeError("Field mismatch: elements belong to different field contexts")

    # ---- core ops  ----
    def _add(self, other: "Element") -> "Element":
        raise NotImplementedError

    def _sub(self, other: "Element") -> "Element":
        raise NotImplementedError

    def _mul(self, other: "Element") -> "Element":
        raise NotImplementedError

    def _neg(self) -> "Element":
        raise NotImplementedError

    def _inv(self) -> "Element":
        raise NotImplementedError

    # ---- operator overloads ----
    def __add__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        self._check_same_field(o)
        return self._add(o)

    def __radd__(self, other: Any) -> "Element":
        return self.__add__(other)

    def __sub__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        self._check_same_field(o)
        return self._sub(o)

    def __rsub__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        o._check_same_field(self)
        return o._sub(self)

    def __mul__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        self._check_same_field(o)
        return self._mul(o)

    def __rmul__(self, other: Any) -> "Element":
        return self.__mul__(other)

    def __truediv__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        self._check_same_field(o)
        return self._mul(o._inv())

    def __rtruediv__(self, other: Any) -> "Element":
        o = self._coerce_other(other)
        o._check_same_field(self)
        return o._mul(self._inv())

    def __neg__(self) -> "Element":
        return self._neg()

    def __pow__(self, e: int) -> "Element":
        if not isinstance(e, int):
            raise TypeError("Exponent must be an int")
        if e < 0:
            return (self._inv()) ** (-e)
        # fast exponentiation
        result = self.field.one
        base = self
        exp = e
        while exp:
            if exp & 1:
                result = result * base
            base = base * base
            exp >>= 1
        return result

    # ---- comparisons / utilities ----
    def is_zero(self) -> bool:
        return self == self.field.zero  # uses __eq__

    def inv(self) -> "Element":
        return self._inv()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Element):
            return False
        if self.field is not other.field:
            return False
        return self.value == other.value

    def __hash__(self) -> int:
        return hash((id(self.field), self.value))
    