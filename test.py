from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

def sortKey(m):
    return m[1][0]

for n in range(5, 6):
    gens = sorted(list(set(Gamma0(n).generators(algorithm="todd-coxeter"))), key = sortKey)
    if n in Primes():
        print(n, len(gens))
        cFs = []
        for gen in gens:
            print(gen, ",")
            cF = contFrac(gen[0][0], gen[1][0])
            #cFs.append(contFrac(gen[0][0], gen[1][0]))
            print(buildMatrix(gen[0][0], gen[1][0]) * T, ",")
            print(cF)
            print("\n")
        #sCFs = sorted(cFs[4:], key = len)
        #for cF in sCFs:
