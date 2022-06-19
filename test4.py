from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

p = 3
G = Gamma0(p)
reps, gens, l, s = G.todd_coxeter()

print(matrix.identity(2), ",")
for i in range(p):
    print(S * T ** i, ",")
print("\n")

for rep in reps:
    print(rep, ",")
print("\n")

for gen in gammaZeroPrimeGens(p):
    print(gen, ",")