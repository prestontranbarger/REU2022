from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

a, c = 17, 2
q1, q2 = 5, 7
dChars1, dChars2 = allDCharacters(q1), allDCharacters(q2)
gamma = buildMatrix(a, c * q1 * q2)
gammaPrime = matrix(ZZ, [[gamma[1][1], -1 * c], [-1 * gamma[0][1] * q1 * q2, a]])
for dChar1 in dChars1:
    for dChar2 in dChars2:
        if (isEven(dChar1) == isEven(dChar2)) and (isPrimitive(dChar1) and isPrimitive(dChar2)):
            lhs, rhs = newFormDedekindSum(dChar1, dChar2, gamma), 0 + 0 * I
            if isEven(dChar1):
                rhs = newFormDedekindSum(dChar2, dChar1, gammaPrime)
            else:
                rhs = -1 * newFormDedekindSum(dChar2, dChar1, gammaPrime) + (1 - psiChar(dChar1, dChar2, gamma)) * LFunction(1, dChar1) * LFunction(1, dChar2) * (gaussSum(conjugateDChar(dChar1)) * gaussSum(conjugateDChar(dChar2)) / (pi * I) ** 2)
            print(float(lhs.real()) + float(lhs.imag()) * I)
            print(float(rhs.real()) + float(rhs.imag()) * I)