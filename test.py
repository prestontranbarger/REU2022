from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

def sortKey(m):
    return m[1][0]

for n in range(3, 25):
    gens = sorted(list(set(Gamma0(n).generators(algorithm="todd-coxeter"))), key = sortKey)
    if n in Primes():
        print(n, len(gens))
        cFs = []
        for gen in gens:
            cFs.append(contFrac(gen[0][0], gen[1][0]))
        sCFs = sorted(cFs[4:], key = len)
        for cF in sCFs:
            print(cF)