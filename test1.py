from SLTwoZ import *
from dedekindSum import *

m = matrix(ZZ, [[17, 32],
                [ 9, 17]])

dChar1 = allDCharacters(3)[1]
dChar2 = allDCharacters(3)[1]
print(isEven(dChar1) == isEven(dChar2))
print(isPrimitive(dChar1))
print(isPrimitive(dChar2))

#for asdf in newFormDedekindSumFast(dChar1, dChar2, m):
#    print(asdf, "\n")

#rwt = TSDecompToRewritingTape(TSDecomp(m))
#print(len(rwt))
#for value in rwt:
#    print(value, "\n")

oneHS = newFormDedekindSum(dChar1, dChar2, m)
print(float(oneHS.real()) + float(oneHS.imag()) * I)
twoHS = newFormDedekindSumFast(dChar1, dChar2, m)
print(float(twoHS.real()) + float(twoHS.imag()) * I)
threeHS = newFormDedekindSumFastPrecompute(dChar1, dChar2, m, chprPathFinder(dChar1, dChar2))
print(float(threeHS.real()) + float(threeHS.imag()) * I)