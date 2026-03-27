from ecc import PrimeField, Curve, group_order, group_structure, point_order

Fp = PrimeField(97)
E = Curve(Fp, a=2, b=3)

P = E.point(3, 6)

print("Curve: y^2 = x^3 + 2x + 3 over F_97")
print()

# --- Point arithmetic ---
print("=== Point Arithmetic ===")
print("P  =", P)
print("2P =", P + P)   # (80,10)
print("3P =", 3 * P)   # (80,87)
print("4P =", 4 * P)   # (3,91) = -P
print("5P =", 5 * P)   # infinity
print()

# --- Group order ---
print("=== Group Order ===")
N = group_order(E)
print(f"#E(F_97) = {N}")
print(f"ord(P)   = {point_order(E, P, N)}")
print()

# --- Group structure ---
print("=== Group Structure ===")
(e1, n1), (e2, n2) = group_structure(E)
print(f"E(F_97) ~= Z/{n1} x Z/{n2}")
print(f"Generator e1 (order {n1}): {e1}")
if n2 > 1:
    print(f"Generator e2 (order {n2}): {e2}")
else:
    print("Group is cyclic (n2 = 1)")