from SLTwoZ import *
from dedekindSum import *
from dirichletCharacters import *

dChar1 = allDCharacters(3)[1]
dChar2 = allDCharacters(5)[2]
q1, q2 = modulus(dChar1), modulus(dChar2)
n = q1 * q2

a, c = 11, 2 * n
gamma = buildMatrix(a, c)

G0G1reps = cosetRepsGammaZeroOverGammaOneGeneralized(n)
inHChecker = inGroupChecker(Gamma1(n))
G0G1 = findCosetRepsLeft(G0G1reps, inHChecker, gamma, n)

print(G0G1 ** (-1) * gamma)
print(TSDecomp(G0G1 ** (-1) * gamma))