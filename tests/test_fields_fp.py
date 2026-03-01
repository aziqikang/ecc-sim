from ecc import PrimeField


def test_fp_basic_arithmetic():
    F = PrimeField(97)

    a = F(5)
    b = F(200)  # coerces mod p -> 6
    assert b.value == 6

    assert a + b == F(11)
    assert a - b == F(96)  # 5 - 6 mod 97
    assert a * b == F(30)
    assert -a + a == F.zero


def test_fp_inverse():
    F = PrimeField(97)
    for v in [1, 2, 3, 5, 11, 42, 96]:
        x = F(v)
        assert x * x.inv() == F.one


def test_fp_inverse_zero_raises():
    F = PrimeField(97)
    try:
        F(0).inv()
        assert False, "Expected ZeroDivisionError"
    except ZeroDivisionError:
        pass
    