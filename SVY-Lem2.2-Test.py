from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

q1, q2 = 5, 7
dChars1, dChars2 = allDCharacters(q1), allDCharacters(q2)
a1, c1 = 11, 5
a2, c2 = 13, 4
gamma1, gamma2 = buildMatrix(a1, c1 * q1 * q2), buildMatrix(a2, c2 * q1 * q2)
for dChar1 in dChars1:
    for dChar2 in dChars2:
        if (isEven(dChar1) == isEven(dChar2)) and (isPrimitive(dChar1) and isPrimitive(dChar2)):
            lhs = newFormDedekindSum(dChar1, dChar2, gamma1 * gamma2)
            rhs = newFormDedekindSum(dChar1, dChar2, gamma1) + psiChar(dChar1, dChar2, gamma1) * newFormDedekindSum(dChar1, dChar2, gamma2)
            print(float(lhs.real()) + float(lhs.imag()) * I)
            print(float(rhs.real()) + float(rhs.imag()) * I)
            print("\n")