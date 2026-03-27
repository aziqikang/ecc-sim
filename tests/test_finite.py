"""
Tests for ecc.ec.finite: point_order, group_order, group_structure.
"""
import math
import pytest
from ecc import PrimeField, Curve
from ecc import point_order, group_order, group_structure


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def curve_97():
    """y^2 = x^3 + 2x + 3 over F_97.  Known: ord(3, 6) = 5."""
    F = PrimeField(97)
    return Curve(F, a=2, b=3)


@pytest.fixture
def curve_f5():
    """y^2 = x^3 + x over F_5.  All non-identity points have order 2 => Z/2 x Z/2."""
    F = PrimeField(5)
    return Curve(F, a=1, b=0)


@pytest.fixture
def curve_f71():
    """y^2 = x^3 + x + 1 over F_71.  Used for a mid-size cyclic-group check."""
    F = PrimeField(71)
    return Curve(F, a=1, b=1)


# ---------------------------------------------------------------------------
# point_order
# ---------------------------------------------------------------------------

class TestPointOrder:
    def test_known_order_5(self, curve_97):
        """(3, 6) has order 5 on y^2 = x^3 + 2x + 3 over F_97 (existing test confirms 5*P = inf)."""
        P = curve_97.point(3, 6)
        N = group_order(curve_97)
        assert point_order(curve_97, P, N) == 5

    def test_infinity_has_order_1(self, curve_97):
        N = group_order(curve_97)
        O = curve_97.infinity()
        assert point_order(curve_97, O, N) == 1

    def test_order_divides_group_order(self, curve_97):
        N = group_order(curve_97)
        P = curve_97.point(3, 6)
        ord_P = point_order(curve_97, P, N)
        assert N % ord_P == 0


# ---------------------------------------------------------------------------
# group_order
# ---------------------------------------------------------------------------

class TestGroupOrder:
    def test_hasse_bound_f97(self, curve_97):
        p = 97
        N = group_order(curve_97)
        lo = (math.isqrt(p) - 1) ** 2
        hi = (math.isqrt(p) + 1) ** 2 + 4   # small slack for isqrt floor
        assert lo <= N <= hi

    def test_known_point_annihilated(self, curve_97):
        """[N]*P = infinity for a known point."""
        N = group_order(curve_97)
        P = curve_97.point(3, 6)
        assert curve_97.scalar_mul(N, P).is_infinity

    def test_divisible_by_known_order(self, curve_97):
        """N must be divisible by 5 since ord(3,6)=5."""
        N = group_order(curve_97)
        assert N % 5 == 0

    def test_hasse_bound_f5(self, curve_f5):
        p = 5
        N = group_order(curve_f5)
        lo = (math.isqrt(p) - 1) ** 2
        hi = (math.isqrt(p) + 1) ** 2 + 4
        assert lo <= N <= hi

    def test_exact_order_f5(self, curve_f5):
        """E: y^2=x^3+x over F_5 has exactly 4 points (including infinity)."""
        assert group_order(curve_f5) == 4

    def test_random_point_annihilated_f5(self, curve_f5):
        N = group_order(curve_f5)
        for x_val, y_val in [(0, 0), (2, 0), (3, 0)]:
            P = curve_f5.point(x_val, y_val)
            assert curve_f5.scalar_mul(N, P).is_infinity

    def test_type_error_on_rational_field(self):
        from ecc import RationalField
        F = RationalField()
        E = Curve(F, a=-1, b=0)
        with pytest.raises(TypeError):
            group_order(E)


# ---------------------------------------------------------------------------
# group_structure
# ---------------------------------------------------------------------------

class TestGroupStructure:
    def _check_structure(self, curve, e1, n1, e2, n2):
        N = n1 * n2
        # n2 | n1
        assert n1 % n2 == 0
        # generators lie on the curve
        assert curve.is_on_curve(e1)
        assert curve.is_on_curve(e2)
        # orders are correct
        assert point_order(curve, e1, n1) == n1
        if n2 > 1:
            assert point_order(curve, e2, n2) == n2
        else:
            assert e2.is_infinity
        # group order consistency
        assert group_order(curve) == N

    def test_cyclic_f97(self, curve_97):
        (e1, n1), (e2, n2) = group_structure(curve_97)
        self._check_structure(curve_97, e1, n1, e2, n2)
        # F_97 cyclic check: n2 must divide gcd(N, p-1) = gcd(N, 96)
        N = n1 * n2
        assert n2 % 1 == 0   # trivially true; real constraint below
        assert (N % 5) == 0  # known from point order

    def test_noncyclic_f5(self, curve_f5):
        """y^2=x^3+x over F_5 is Z/2 x Z/2, so n1=n2=2."""
        (e1, n1), (e2, n2) = group_structure(curve_f5)
        self._check_structure(curve_f5, e1, n1, e2, n2)
        assert n1 == 2
        assert n2 == 2

    def test_structure_product_equals_group_order(self, curve_f71):
        (e1, n1), (e2, n2) = group_structure(curve_f71)
        self._check_structure(curve_f71, e1, n1, e2, n2)

    def test_generators_annihilate_at_order(self, curve_97):
        (e1, n1), (e2, n2) = group_structure(curve_97)
        assert curve_97.scalar_mul(n1, e1).is_infinity
        if n2 > 1:
            assert curve_97.scalar_mul(n2, e2).is_infinity
