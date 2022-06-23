from additiveCharacters import *
from dirichletCharacters import *
from dedekindSum import *
from SLTwoZ import *

q1, q2 = 9, 27
a, c = 23, q1 * q2
dChar1 = allDCharacters(q1)[1]
dChar2 = allDCharacters(q2)[5]
gamma = buildMatrix(a, c)

print(isEven(dChar1) == isEven(dChar2))
print(isPrimitive(dChar1))
print(isPrimitive(dChar2))

nFDS = newFormDedekindSum(dChar1, dChar2, gamma)
print(float(nFDS.real()) + float(nFDS.imag()) * I)

nFDs = newFormDedekindSumNaiveFast(dChar1, dChar2, gamma)
print(float(nFDS.real()) + float(nFDS.imag()) * I)