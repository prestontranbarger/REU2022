from dedekindSum import *

dChar1 = allDCharacters(5)[3]
dChar2 = allDCharacters(7)[5]
print(isEven(dChar1))
print(isEven(dChar1) == isEven(dChar2))
print(isPrimitive(dChar1))
print(isPrimitive(dChar2))

m = matrix(ZZ, [[351, 269],
                [595, 456]])

print(inGammaOne(m, modulus(dChar1) * modulus(dChar2)))

lhs = newFormDedekindSum(dChar1, dChar2, m)
rhs = newFormDedekindSumSVYThm1Dot2(dChar1, dChar2, m)

print(lhs.real() + lhs.imag() * I)
print(rhs.real() + rhs.imag() * I)