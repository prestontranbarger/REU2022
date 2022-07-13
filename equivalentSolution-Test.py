from additiveCharacters import *
from dirichletCharacters import *
from dedekindSum import *
from SLTwoZ import *

q1, q2 = 3, 15
a, c = 29, 4 * q1 * q2
dChar1 = allDCharacters(q1)[1]
dChar2 = allDCharacters(q2)[5]
gamma = buildMatrix(a, c)

print(isEven(dChar1) == isEven(dChar2))
print(isPrimitive(dChar1))
#for i in range(q1):
#    print(i, dChar1(i))
print(isPrimitive(dChar2))
#for j in range(q2):
#    print(j, dChar2(j))
chprPath = chprPathFinder(dChar1, dChar2)
print(chprPath)

nFDS = newFormDedekindSum(dChar1, dChar2, gamma)
print(float(nFDS.real()) + float(nFDS.imag()) * I)

#nFDs = newFormDedekindSumNaiveFast(dChar1, dChar2, gamma)
#print(float(nFDS.real()) + float(nFDS.imag()) * I)

#nFDs = newFormDedekindSumFast(dChar1, dChar2, gamma)
#print(float(nFDS.real()) + float(nFDs.imag()) * I)

nFDs = newFormDedekindSumNaiveFastGeneralizedPrecompute(dChar1, dChar2, gamma, chprPath[0])
print(float(nFDS.real()) + float(nFDs.imag()) * I)

nFDs = newFormDedekindSumFastGeneralizedPrecompute(dChar1, dChar2, gamma, chprPath[1])
print(float(nFDS.real()) + float(nFDs.imag()) * I)