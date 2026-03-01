from ecc import PrimeField, Curve


def test_associativity_spotcheck_fp97():
    Fp = PrimeField(97)
    E = Curve(Fp, a=2, b=3)

    P = E.point(3, 6)
    Q = 2 * P
    R = 3 * P

    left = (P + Q) + R
    right = P + (Q + R)
    assert left == right
