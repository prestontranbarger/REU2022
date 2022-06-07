from dirichletCharacters import *

def j(gamma, z):
    return gamma[1][0] * z + gamma[1][1]

def gammaOfZ(gamma, z):
    return (gamma[0][0] * z + gamma[0][1]) / (gamma[1][0] * z + gamma[1][1])

#def dedekindEta(z, t = 100):
#    q = e ** (2 * pi * I * z)
#    p = qexp_eta(ZZ[['x']], t).polynomial()
#    print(p)
#    return q ** (1 / 24) * p(x = q)

def dedekindEpsilon(gamma):
    #computes the dedekind epsilon function of a, b, c, d
    a, b, c, d = gamma[0][0], gamma[0][1], gamma[1][0], gamma[1][1]
    if c == 0 and d == 1:
        return e ** (pi * I * b / 12)
    return e ** (pi * I * ((a + d) / (12 * c) - dedekindSum(d, c) - 1 / 4))

def partition(n):
    #computes the number of partitions of [0, 1, ..., n], returns [p(0), p(1), ..., p(n)]
    parts = [1] + [0] * n
    for t in range(1, n + 1):
        for i, x in enumerate(range(t, n + 1)):
            parts[x] += parts[i]
    return parts

def sawtooth(x):
    #sawtooth function
    if x == floor(x):
        return 0
    return x - floor(x) - 1 / 2

def dedekindSum(b, c):
    sum = 0
    i = 0
    while b != 0 and c != 1:
        b %= c
        sum += (-1) ** i * ((b ** 2 + c ** 2 + 1) / (12 * b * c) - 1 / 4)
        b, c = c, b
        i += 1
    return sum

def dedekindSumSlow(b, c):
    return generalizedDedekindSum(1, b, c)

def generalizedDedekindSum(a, b, c):
    sum = 0
    for n in range(c):
        sum += sawtooth(a * n / c) * sawtooth(b * n / c)
    return sum

def newFormDedekindSum(dChar1, dChar2, gamma):
    sum = 0
    q1 = modulus(dChar1)
    a, c = gamma[0][0], gamma[1][0]
    dictChar1 = getValuesDict(dChar1)
    dictChar2 = getValuesDict(dChar2)
    for j in range(c):
        for n in range(modulus(dChar1)):
            sum += dictChar2(j).conjugate() * dictChar1(n).conjugate() * sawtooth(j / c) * sawtooth(n / q1 + a * j / c)
    return sum