from dirichletCharacters import *
from sage.libs.lcalc.lcalc_Lfunction import *

def j(gamma, z):
    return gamma[1][0] * z + gamma[1][1]

def gammaOfZ(gamma, z):
    return (gamma[0][0] * z + gamma[0][1]) / (gamma[1][0] * z + gamma[1][1])

def dedekindEta(z, t = 100):
    #TODO: needs fixing
    q = e ** (2 * pi * I * z)
    p = qexp_eta(ZZ[['x']], t).polynomial()
    print(p)
    return q ** (1 / 24) * p(x = q)

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

def psiChar(dChar1, dChar2, gamma):
    d = gamma[1][1]
    return dChar1(d) * dChar2(d).conjugate()

def LFunction(s, dChar):
    return Lfunction_from_character(dChar).value(s)

def dedekindSum(b, c):
    #computes the standard dedekind sum using reciprocity laws to speed up computation
    #see: https://en.wikipedia.org/wiki/Dedekind_sum
    sum = 0
    i = 0
    while b != 0 and c != 1:
        b %= c
        sum += (-1) ** i * ((b ** 2 + c ** 2 + 1) / (12 * b * c) - 1 / 4)
        b, c = c, b
        i += 1
    return sum

def dedekindSumSlow(b, c):
    # computes the standard dedekind sum
    # see: https://en.wikipedia.org/wiki/Dedekind_sum
    return generalizedDedekindSum(1, b, c)

def generalizedDedekindSum(a, b, c):
    # computes the generalized dedekind sum
    # see: https://en.wikipedia.org/wiki/Dedekind_sum
    sum = 0
    for n in range(c):
        sum += sawtooth(a * n / c) * sawtooth(b * n / c)
    return sum

def newFormDedekindSum(dChar1, dChar2, gamma):
    sum = 0
    q1, q2 = modulus(dChar1), modulus(dChar2)
    a, c = gamma[0][0], gamma[1][0]
    #print(gamma)
    for j in range(abs(c)):
        for n in range(q1):
            #print(sum)
            sum += dChar2(j).conjugate() * dChar1(n).conjugate() * sawtooth(j / c) * sawtooth(n / q1 + a * j / c)
            #print(sum)
    #print(sum)
    return sum

def cdMaxNorm(norm = 10000):
    pairs = []
    for c in range(ceil(-1 * sqrt(norm)), floor(sqrt(norm))):
        for d in range(ceil(-1 * sqrt(norm - c ** 2)), floor(sqrt(norm - c ** 2))):
            if gcd(c, d) == 1:
                pairs.append([c, d])
    return pairs

def newformEisensteinSeries(dChar1, dChar2, z, s, norm = 10000):
    #TODO: needs fixing
    if isEven(dChar1) == isEven(dChar2):
        q1, q2 = modulus(dChar1), modulus(dChar2)
        sum = 0
        for pair in cdMaxNorm(norm):
            c, d = pair[0], pair[1]
            sum += ((q2 * z.imag()) ** s * dChar1(c) * dChar2(d)) / (c * q2 * z + d) ** (2 * s)
        return sum / 2
    return None

def completedEisensteinSeries(dChar1, dChar2, z, s, norm = 10000):
    #TODO: needs fixing
    if isEven(dChar1) == isEven(dChar2):
        q1, q2 = modulus(dChar1), modulus(dChar2)
        return ((q2 / pi) ** s / gaussSum(dChar2)) * gammaFunction(s) * LFunction(2 * s, prodDCharacters(dChar1, dChar2)) * newformEisensteinSeries(dChar1, dChar2, z, s, norm)
    return None