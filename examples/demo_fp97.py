from ecc import PrimeField, Curve

Fp = PrimeField(97)
E = Curve(Fp, a=2, b=3)

P = E.point(3, 6)

print("Curve: y^2 = x^3 + 2x + 3 over F_97")
print("P =", P)
print("2P =", P + P)   # (80,10)
print("3P =", 3 * P)   # (80,87)
print("4P =", 4 * P)   # (3,91) = -P
print("5P =", 5 * P)   # infinity