# ecc/ec/finite.py
"""
Algorithms for elliptic curves over finite (prime) fields.

Public API:
    point_order(curve, P, N)     -> int
    group_order(curve)           -> int
    group_structure(curve)       -> ((Point, int), (Point, int))
"""
from __future__ import annotations

import math
import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .curve import Curve
    from .point import Point


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_prime_field(curve: "Curve") -> int:
    """Return p if curve is over PrimeField(p), else raise TypeError."""
    p = getattr(curve.field, "p", None)
    if not isinstance(p, int):
        raise TypeError(
            "group_order and group_structure require a curve over a PrimeField"
        )
    return p


def _factor(n: int) -> dict[int, int]:
    """Prime factorization via trial division. Returns {prime: exponent}."""
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = 1
    return factors


def _tonelli_shanks(n: int, p: int) -> int:
    """Modular square root: return y with y^2 == n (mod p). n must be a QR mod p."""
    if n % p == 0:
        return 0
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    # General Tonelli-Shanks
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    M, c, t, R = S, pow(z, Q, p), pow(n, Q, p), pow(n, (Q + 1) // 2, p)
    while True:
        if t == 1:
            return R
        i, tmp = 1, t * t % p
        while tmp != 1:
            tmp = tmp * tmp % p
            i += 1
        b = pow(c, 1 << (M - i - 1), p)
        M, c, t, R = i, b * b % p, t * b * b % p, R * b % p


def _random_point(curve: "Curve", p: int) -> "Point":
    """Return a uniformly random affine point on curve over F_p."""
    a_val = curve.a.value
    b_val = curve.b.value
    while True:
        x_val = random.randrange(p)
        rhs = (pow(x_val, 3, p) + a_val * x_val + b_val) % p
        if rhs == 0:
            return curve.point(x_val, 0)
        if pow(rhs, (p - 1) // 2, p) != 1:
            continue  # not a quadratic residue
        y_val = _tonelli_shanks(rhs, p)
        if random.randrange(2):
            y_val = p - y_val
        return curve.point(x_val, y_val)


def _p_val(n: int, p: int) -> int:
    """p-adic valuation v_p(n)."""
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


# ---------------------------------------------------------------------------
# point_order
# ---------------------------------------------------------------------------

def point_order(curve: "Curve", P: "Point", N: int) -> int:
    """
    Return ord(P), given that ord(P) divides N.

    Uses the divide-out algorithm: for each prime power in the factorization
    of N, repeatedly divide out the prime as long as the result still
    annihilates P. Complexity: O(log^2 N) scalar multiplications.
    """
    factors = _factor(N)
    order = N
    for prime, exp in factors.items():
        for _ in range(exp):
            candidate = order // prime
            if curve.scalar_mul(candidate, P).is_infinity:
                order = candidate
            else:
                break
    return order


# ---------------------------------------------------------------------------
# group_order
# ---------------------------------------------------------------------------

# Below this threshold, direct point counting is used instead of BSGS.
# For small p the group exponent may be <= the BSGS step size m, causing
# all random points to have small order and the baby-step table to produce
# wrong indices.  Direct counting is O(p) and always exact.
_DIRECT_COUNT_THRESHOLD = 1_000


def group_order(curve: "Curve") -> int:
    """
    Compute #E(F_p).

    For p < 1000 uses direct point counting (O(p), exact).
    For p >= 1000 uses Shanks' baby-step giant-step (O(p^(1/4))).

    Hasse's theorem bounds N = #E(F_p) to the interval
        (sqrt(p) - 1)^2  <=  N  <=  (sqrt(p) + 1)^2.

    Requires curve.field to be a PrimeField.
    """
    p = _require_prime_field(curve)
    if p < _DIRECT_COUNT_THRESHOLD:
        return _count_points(curve, p)
    return _bsgs_group_order(curve, p)


def _count_points(curve: "Curve", p: int) -> int:
    """Count #E(F_p) by iterating over all x in F_p. O(p), exact."""
    a_val = curve.a.value
    b_val = curve.b.value
    N = 1  # point at infinity
    for x in range(p):
        rhs = (pow(x, 3, p) + a_val * x + b_val) % p
        if rhs == 0:
            N += 1                    # one affine point (x, 0)
        elif pow(rhs, (p - 1) // 2, p) == 1:
            N += 2                    # two affine points (x, ±y)
    return N


def _bsgs_group_order(curve: "Curve", p: int) -> int:
    """
    Baby-step giant-step for p >= _DIRECT_COUNT_THRESHOLD.

    Write N = p+1-t, |t| <= B = ceil(2*sqrt(p)), shift to u = t+B in [0,2B],
    and factor u = j + k*m (m ~ sqrt(2B)).  Collision [p+1+B-km]*P = [j]*P
    gives candidate N = p+1-(u-B), verified with an independent random point.

    If P has order <= m, the baby-step table has overwritten entries yielding
    wrong j values.  This is detected during table construction and the point
    is discarded.
    """
    B = 2 * math.isqrt(p) + 2
    m = math.isqrt(2 * B) + 2
    lo = max(1, p + 1 - B)
    hi = p + 1 + B

    while True:
        P = _random_point(curve, p)

        # Baby steps: detect small-order P (order <= m) by collision in table.
        baby: dict = {}
        acc = curve.infinity()
        small_order = False
        for j in range(m + 1):
            if acc in baby:
                small_order = True
                break
            baby[acc] = j
            acc = curve.add(acc, P)
        if small_order:
            continue

        M = curve.scalar_mul(m, P)
        neg_M = curve.neg(M)
        Q_shift = curve.scalar_mul(p + 1 + B, P)

        # Giant steps: R_k = [(p+1+B) - k*m]*P
        candidates: list[int] = []
        R = Q_shift
        for k in range(2 * m + 3):
            if R in baby:
                u = baby[R] + k * m
                t = u - B
                N_cand = p + 1 - t
                if lo <= N_cand <= hi:
                    candidates.append(N_cand)
            R = curve.add(R, neg_M)

        # Verify each candidate with an independent random point
        for N_cand in dict.fromkeys(candidates):
            Q = _random_point(curve, p)
            if curve.scalar_mul(N_cand, Q).is_infinity:
                return N_cand
        # All candidates failed verification; retry with a new point.


# ---------------------------------------------------------------------------
# group_structure
# ---------------------------------------------------------------------------

def group_structure(
    curve: "Curve",
) -> tuple[tuple["Point", int], tuple["Point", int]]:
    """
    Compute the group structure of E(F_p).

    By the structure theorem for finite abelian groups,
        E(F_p) ~= Z/n1 x Z/n2,
    where n2 | n1 and n1 * n2 = N = #E(F_p).  The Weil pairing forces
    n2 | gcd(N, p-1), which is a strong constraint for large prime fields.

    Returns ((e1, n1), (e2, n2)) where:
      - e1 has order n1 and e2 has order n2
      - e1 and e2 generate E(F_p)
      - if the group is cyclic (n2 = 1), e2 is the point at infinity

    Algorithm:
      1. Compute N via baby-step giant-step.
      2. For each prime q with q^a || N, find r = max v_q(ord(Q)) over random
         q-primary points Q.  Set s = a - r; then the q-primary component is
         Z/q^r x Z/q^s.
      3. Combine across primes: n1 = prod q^r, n2 = prod q^s.
      4. Find generators e1 (order n1) and e2 (order n2, independent of e1).

    Requires curve.field to be a PrimeField.
    """
    p = _require_prime_field(curve)
    N = group_order(curve)
    factors = _factor(N)

    n1 = 1
    n2 = 1
    for q, a in factors.items():
        r = _max_q_exponent(curve, N, q, a, p)
        s = a - r
        n1 *= q ** r
        n2 *= q ** s

    e1 = _find_generator(curve, N, n1, p)

    if n2 == 1:
        e2 = curve.infinity()
    else:
        e2 = _find_second_generator(curve, N, n1, n2, e1, factors, p)

    return (e1, n1), (e2, n2)


# ---------------------------------------------------------------------------
# group_structure internals
# ---------------------------------------------------------------------------

def _max_q_exponent(
    curve: "Curve", N: int, q: int, a: int, p: int
) -> int:
    """
    Find r = max { v_q(ord(Q)) : Q a q-primary point in E(F_p) }.

    With high probability over random sampling this equals the exponent r of
    the q-primary component Z/q^r x Z/q^s (r >= s).  Once r is found,
    s = a - r follows immediately.
    """
    cofactor = N // (q ** a)
    max_r = 0
    for _ in range(max(24, 6 * a)):
        P = _random_point(curve, p)
        Q = curve.scalar_mul(cofactor, P)
        if Q.is_infinity:
            continue
        # Divide out q until [order/q]*Q != O: gives exact q-power order of Q.
        order = q ** a
        for _ in range(a):
            if curve.scalar_mul(order // q, Q).is_infinity:
                order //= q
            else:
                break
        r = _p_val(order, q)
        if r > max_r:
            max_r = r
            if max_r == a:      # can't improve further
                break
    return max_r


def _find_generator(curve: "Curve", N: int, n: int, p: int) -> "Point":
    """
    Return a point of order exactly n on curve (n | N).

    Samples random points directly and checks their order against n.
    Using cofactor = N//n to pre-filter is wrong when gcd(cofactor, n) > 1:
    the image [cofactor]*E is a proper subgroup and may not contain
    any element of order n (e.g. for E ~= Z/50 x Z/2, cofactor=2 maps
    into Z/25, which has no element of order 50 = n1).
    """
    while True:
        P = _random_point(curve, p)
        if point_order(curve, P, N) == n:
            return P


def _find_second_generator(
    curve: "Curve",
    N: int,
    n1: int,
    n2: int,
    e1: "Point",
    factors: dict[int, int],
    p: int,
) -> "Point":
    """
    Find e2 with ord(e2) = n2, independent from e1.

    Works prime-by-prime over the factorization of n2, then combines via CRT:

      For each prime q with q^s || n2 and q^r || n1  (r + s = v_q(N) = a, r >= s):
        - cofactor_q = N / q^a maps a random point into the q^a-torsion.
        - g1_q = [n1/q^r]*e1 has order q^r and generates the first cyclic
          q-primary factor Z/q^r.
        - Iterate k = 0..q^r-1: try acc = Q - k*g1_q.  Among these, exactly
          the k-values with alpha(Q) - k ≡ 0 (mod q^{r-s}) (where alpha(Q)
          is Q's first-component DL) yield ord(acc) = q^s AND acc outside
          <g1_q>.  Return the first such acc as gen_q.

      CRT recombination: e2 = sum_q  [n2/q^s] * gen_q.
      Each term has order q^s; since the primes are distinct the terms are
      in orthogonal subgroups and e2 has order lcm(q^s) = n2.

    Uses: cofactor = N // q^a  (NOT N // n2, which equals n1 and kills all points
    because the group exponent is n1).
    """
    e2 = curve.infinity()
    for q, a in factors.items():
        s = _p_val(n2, q)
        if s == 0:
            continue
        gen_q = _find_q_primary_gen2(curve, N, n1, n2, e1, q, a, s, p)
        crt_coeff = n2 // (q ** s)              # kills all other q'-primary parts
        e2 = curve.add(e2, curve.scalar_mul(crt_coeff, gen_q))
    return e2


def _find_q_primary_gen2(
    curve: "Curve",
    N: int, n1: int, n2: int, e1: "Point",
    q: int, a: int, s: int, p: int,
) -> "Point":
    """
    Find a point of order q^s in the second cyclic q-primary factor.

    g1_q = [n1/q^r]*e1  (order q^r, generates the first cyclic q-factor).
    For random Q = [N/q^a]*P, we try acc = Q - k*g1_q for k = 0..q^r-1.
    When acc is outside <g1_q> and has order q^s, it generates the second
    cyclic q-factor.  Expected iterations: O(q^r) per P, O(1) P samples.
    """
    r = a - s
    cofactor_q = N // (q ** a)
    g1_q = curve.scalar_mul(n1 // (q ** r), e1)   # order q^r
    neg_g1_q = curve.neg(g1_q)
    target_order = q ** s
    q_r = q ** r

    while True:
        P = _random_point(curve, p)
        Q = curve.scalar_mul(cofactor_q, P)
        if Q.is_infinity:
            continue
        acc = Q
        for _ in range(q_r):
            if not acc.is_infinity:
                # Check order == q^s
                ord_acc = q ** a
                for _ in range(a):
                    if curve.scalar_mul(ord_acc // q, acc).is_infinity:
                        ord_acc //= q
                    else:
                        break
                # Check independence: acc not in <g1_q>
                if ord_acc == target_order and not _in_cyclic_span(curve, acc, g1_q, q_r):
                    return acc
            acc = curve.add(acc, neg_g1_q)   # acc = Q - (k+1)*g1_q


def _in_cyclic_span(
    curve: "Curve", point: "Point", generator: "Point", order: int
) -> bool:
    """Return True if point is in <generator>, where generator has given order."""
    acc = curve.infinity()
    for _ in range(order):
        if acc == point:
            return True
        acc = curve.add(acc, generator)
    return False
