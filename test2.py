from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

def sortKey(m):
    return m[1][0]

p = 5
gens = sorted(list(set(Gamma0(p).generators(algorithm="todd-coxeter"))), key = sortKey)
for gen in gens:
    print(gen, ",")

print("\n")

for gen in gammaZeroPrimeGens(p):
    print(gen, ",")