# ECC

Toy elliptic curve playground with field-agnostic arithmetic.

## Currently Supported Fields
- $`\mathbb{Q}`$ (rationals via `fractions.Fraction`)
- $`\mathbb{F}_p`$ (prime fields)

## Curve Model
Short Weierstrass form
```math
y^2 = x^3 + ax + b
```


> Characteristic 2 and 3 are not supported.

---

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate

pip install -r requirements-dev.txt
pip install -e .
```

---

## Run Tests

```bash
python -m pytest -q
```

---

## Run Examples

```bash
python -m examples.demo_fp97
python -m examples.demo_q
python -m examples.demo_group_structure
```

---

## Project Structure

```
ecc/        library code
tests/      pytest tests
examples/   demo scripts
```

---

## Mathematical Background

### Elliptic curves over finite fields

For a prime $p > 3$, the group $E(\mathbb{F}_p)$ of points on a short Weierstrass curve $y^2 = x^3 + ax + b$ over $\mathbb{F}_p$ is a finite abelian group under the chord-and-tangent law.  By the structure theorem for finite abelian groups:

```math
E(\mathbb{F}_p) \;\cong\; \mathbb{Z}/n_1\mathbb{Z} \;\times\; \mathbb{Z}/n_2\mathbb{Z}
```

where $n_2 \mid n_1$ and $n_1 n_2 = N = \#E(\mathbb{F}_p)$.  The Weil pairing forces $n_2 \mid \gcd(N,\, p-1)$, a strong constraint that often makes the group cyclic ($n_2 = 1$) for random curves over large prime fields.

The representation used in this library is $((e_1,\, n_1),\, (e_2,\, n_2))$, where $e_1$ and $e_2$ are explicit points of orders $n_1$ and $n_2$ that generate the two cyclic factors.

---

### Step 1 — Computing $N = \#E(\mathbb{F}_p)$

**Hasse's theorem** gives the key bound:

```math
\bigl|\, N - (p+1) \,\bigr| \;\leq\; 2\sqrt{p}
\qquad\Longrightarrow\qquad
N = p + 1 - t, \quad |t| \leq 2\sqrt{p}
```

The integer $t$ is the *trace of Frobenius*.  Hasse's bound confines $N$ to a window of width $\approx 4\sqrt{p}$ centered at $p+1$, which the BSGS algorithm exploits.

**For small $p$ (< 1000):** points are counted directly in $O(p)$ by iterating $x \in \mathbb{F}_p$ and checking whether $x^3 + ax + b$ is a quadratic residue.

**For large $p$:** Shanks' *baby-step giant-step* (BSGS) runs in $O(p^{1/4})$ time and space:

1. Fix $B = \lceil 2\sqrt{p}\,\rceil$ and step size $m = \lceil\sqrt{2B}\,\rceil$ so that $m^2 > 2B$.
2. Shift the unknown $t$ to $u = t + B \in [0,\,2B]$ and write $u = j + km$ with $j \in [0,m)$ and $k \in [0, 2m]$.
3. **Baby steps.** Build the table $\{[j]P \mapsto j\}$ for $j = 0, \ldots, m$.  If a collision appears before $j = m$, $P$ has order $\leq m$ and is discarded.
4. **Giant steps.** Starting from $R_0 = [(p+1+B)]P$, iterate $R_k = R_{k-1} - [m]P$.  A match $R_k = [j]P$ gives $u = j + km$, hence $t = u - B$ and the candidate $N = p + 1 - t$.
5. Verify each candidate $N_\text{cand}$ with an independent random point $Q$: accept if $[N_\text{cand}]Q = \mathcal{O}$.

---

### Step 2 — Finding the group structure $(n_1, n_2)$

Given $N$ (and its factorization $N = \prod q_i^{a_i}$), the decomposition proceeds **prime by prime**.

For each prime $q$ with $q^a \| N$, the $q$-primary component is $\mathbb{Z}/q^r \times \mathbb{Z}/q^s$ with $r \geq s \geq 0$ and $r + s = a$.  To find $r$:

- Sample random points and multiply by the cofactor $N/q^a$ to land in the $q^a$-torsion.
- Find each point's exact $q$-power order by repeatedly dividing out $q$.
- The maximum $q$-power order seen over sufficiently many samples converges to $q^r$ with high probability.

Then $s = a - r$, and across all primes:
```math
n_1 = \prod q^r, \qquad n_2 = \prod q^s.
```

**Finding $e_1$ (order $n_1$).** Sample random $P \in E(\mathbb{F}_p)$ until $\mathrm{ord}(P) = n_1$.

**Finding $e_2$ (order $n_2$, independent of $e_1$).** For each prime $q \mid n_2$ (with $q^s \| n_2$, $q^r \| n_1$):
1. Let $g_{1,q} = [n_1/q^r]\,e_1$ (order $q^r$, generates the first cyclic $q$-factor).
2. Pick random $P$, compute $Q = [N/q^a]\,P$ (in the $q^a$-torsion).
3. Iterate $\mathrm{acc} = Q - k\,g_{1,q}$ for $k = 0, \ldots, q^r - 1$; return the first $\mathrm{acc}$ with $\mathrm{ord}(\mathrm{acc}) = q^s$ and $\mathrm{acc} \notin \langle g_{1,q} \rangle$.

> **Why not use cofactor $N/n_2$?**  Because $N/n_2 = n_1$ is the group exponent, so $[n_1]P = \mathcal{O}$ for *every* $P$ — the image is trivially $\{\mathcal{O}\}$ and no non-trivial second generator can be found.  The $q$-primary cofactor $N/q^a$ correctly targets the joint torsion without over-killing.

The per-prime generators are combined into a single $e_2$ of order $n_2$ via the **Chinese Remainder Theorem**:
```math
e_2 \;=\; \sum_{q \mid n_2} \Bigl[\frac{n_2}{q^s}\Bigr]\, \mathrm{gen}_q
```
where each term has order $q^s$ and the terms are supported on disjoint primes, giving $\mathrm{ord}(e_2) = \mathrm{lcm}(q^s) = n_2$.
