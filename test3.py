from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

p = 13
dChar1 = allDCharacters(1)[0]
dChar2 = allDCharacters(p)[2]
abVals = []
if (isEven(dChar1) == isEven(dChar2)) and (isPrimitive(dChar1) and isPrimitive(dChar2)):
    for gen in gammaZeroPrimeGens(p):
        #print(contFrac(gen[0][0], gen[1][0]))
        nFDS = newFormDedekindSum(dChar1, dChar2, gen)
        ab = abs(nFDS)
        print(psiChar(dChar1, dChar2, gen) * nFDS)
        #print(nFDS, ab)
        abVals.append(ab)
    #print(sorted(abVals)[3:])
else:
    print("check dChars")