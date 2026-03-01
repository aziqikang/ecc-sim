# ECC

Toy elliptic curve playground with field-agnostic arithmetic.

## Currently Supported Fields
- $`\mathbb{Q}`$ (rationals via `fractions.Fraction`)
- $`\mathbb{F}_p`$ (prime fields)

## Curve Model
Short Weierstrass form
$$ y^2 = x^3 + ax + b $$

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
```

---

## Project Structure

```
ecc/        library code
tests/      pytest tests
examples/   demo scripts
```