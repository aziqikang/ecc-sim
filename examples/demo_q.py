from ecc import RationalField, Curve

Q = RationalField()
E = Curve(Q, a=1, b=1)  # y^2 = x^3 + x + 1

P = E.point(0, 1)
O = E.infinity()

print("Curve: y^2 = x^3 + x + 1 over Q")
print("P =", P)
print("O =", O)
print("2P =", P + P)
print("3P =", 3 * P)
print("P + (-P) =", P + (-P))
