from ecc import RationalField, Curve


def test_curve_over_q_identity_inverse_and_closure():
    Q = RationalField()
    E = Curve(Q, a=1, b=1)  # y^2 = x^3 + x + 1
    O = E.infinity()

    P = E.point(0, 1)
    assert E.is_on_curve(P)

    assert P + O == P
    assert O + P == P
    assert (P + (-P)).is_infinity

    S = P + P
    assert E.is_on_curve(S)
    assert not S.is_infinity


def test_curve_over_q_y_zero_doubling_is_infinity():
    # y^2 = x^3 - x has points (0,0), (1,0), (-1,0); these are 2-torsion so doubling -> infinity
    Q = RationalField()
    E = Curve(Q, a=-1, b=0)
    O = E.infinity()

    P = E.point(0, 0)
    assert (P + P) == O
    