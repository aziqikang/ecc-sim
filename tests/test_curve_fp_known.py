from ecc import PrimeField, Curve


def test_curve_over_fp_known_values_mod_97():
    Fp = PrimeField(97)
    E = Curve(Fp, a=2, b=3)

    P = E.point(3, 6)
    assert E.is_on_curve(P)

    twoP = P + P
    assert repr(twoP) == "Point(Fp_97(80), Fp_97(10))"

    threeP = twoP + P
    assert repr(threeP) == "Point(Fp_97(80), Fp_97(87))"

    fourP = 4 * P
    assert repr(fourP) == "Point(Fp_97(3), Fp_97(91))"

    fiveP = 5 * P
    assert fiveP.is_infinity

    O = E.infinity()
    assert P + O == P
    assert O + P == P
    assert (P + (-P)).is_infinity
