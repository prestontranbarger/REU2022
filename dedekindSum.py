from base import *

def j(c, d, z):
    return c * z + d

def gamma(a, b, c, d, z):
    return (a * z + b) / (c * z + d)

#def dedekindEta(z):

def dedekindEpsilon(a, b, c, d):
    #computes the dedekind epsilon function of a, b, c, d
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