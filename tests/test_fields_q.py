from fractions import Fraction
from ecc import RationalField


def test_q_coercion_and_repr():
    Q = RationalField()
    a = Q(2) / Q(3)
    assert repr(a) == "Q(2/3)"
    assert Q(Fraction(2, 3)) == a


def test_q_field_axioms_spot_checks():
    Q = RationalField()
    a = Q(2) / Q(3)
    b = Q(5) / Q(7)

    assert a + Q.zero == a
    assert a * Q.one == a
    assert a - a == Q.zero
    assert (-a) + a == Q.zero

    assert (a + b) - b == a
    assert (a * b) / b == a

    assert (a**0) == Q.one
    assert (a**3) == a * a * a


def test_q_inverse_zero_raises():
    Q = RationalField()
    try:
        Q.zero.inv()
        assert False, "Expected ZeroDivisionError"
    except ZeroDivisionError:
        pass
