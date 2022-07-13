from dirichletCharacters import *
from SLTwoZ import *
from precomputation import *

from sage.libs.lcalc.lcalc_Lfunction import *

def j(gamma, z):
    #given a matrix gamma in SL2Z, this returns j(z)=c*z+d
    return gamma[1][0] * z + gamma[1][1]

def gammaOfZ(gamma, z):
    #given a matrix gamma in SL2Z, this returns the LFT gamma*z=(az+b)/(cz+d)
    return (gamma[0][0] * z + gamma[0][1]) / (gamma[1][0] * z + gamma[1][1])

#def dedekindEta(z, t = 100):
#    #TODO: needs fixing
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

def psiChar(dChar1, dChar2, gamma):
    #given two characters and a matrix dChar1, dChar2, and gamma respectively, this returns dChar1(d)*conjugate(dChar2(d))
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
    #computes the new form dedekind sum of gamma given two primative characters with similar parity
    #this is in accordance with SVY's definition of a finite double sum formula
    sum = 0
    q1, q2 = modulus(dChar1), modulus(dChar2)
    a, c = gamma[0][0] if gamma[1][0] > 0 else -1 * gamma[0][0],\
           gamma[1][0] if gamma[1][0] > 0 else -1 * gamma[1][0]
    for j in range(c):
        for n in range(q1):
            sum += dChar2(j).conjugate() * dChar1(n).conjugate() * (sawtooth(j / c) + 0 * I) * (sawtooth(n / q1 + a * j / c) + 0 * I)
    return sum

def newFormDedekindSumSVYThm1Dot2(dChar1, dChar2, gamma):
    # computes the new form dedekind sum of gamma given two primative characters with similar parity
    # this is in accordance with SVY's definition of a finite double sum formula
    sum = 0
    q1, q2 = modulus(dChar1), modulus(dChar2)
    if abs(gamma[0][1] * q2 ** 2 * q1) < abs(gamma[1][0] * q1):
        at, ct = gamma[1][1], -1 * gamma[0][1] * q1 * q2
        a, c = at if ct > 0 else -1 * at,\
               ct if ct > 0 else -1 * ct
        for j in tqdm(range(c)):
            for n in range(q2):
                sum += dChar1(j).conjugate() * dChar2(n).conjugate() * (sawtooth(j / c) + 0 * I) * (sawtooth(n / q2 + a * j / c) + 0 * I)
        if not isEven(dChar1):
            sum = -1 * sum + (1 - psiChar(dChar1, dChar2, gamma)) * LFunction(1, dChar1) * LFunction(1, dChar2) * gaussSum(conjugateDChar(dChar1)) * gaussSum(conjugateDChar(dChar2)) / ((pi * I) ** 2)
    else:
        a, c = gamma[0][0] if gamma[1][0] > 0 else -1 * gamma[0][0], \
               gamma[1][0] if gamma[1][0] > 0 else -1 * gamma[1][0]
        for j in tqdm(range(c)):
            for n in range(q1):
                sum += dChar2(j).conjugate() * dChar1(n).conjugate() * (sawtooth(j / c) + 0 * I) * (sawtooth(n / q1 + a * j / c) + 0 * I)
    return sum

def newFormDedekindSumNaiveFast(dChar1, dChar2, gamma):
    # computes the new form dedekind sum of gamma given two primative characters with similar parity
    # uses naive fast technique developed during 22 REU
    q1, q2 = modulus(dChar1), modulus(dChar2)
    pppc = pOppOc(q1 * q2)
    #if pppc[0] == 'p':
    #    return newFormDedekindSumNaiveFastPrime(dChar1, dChar2, gamma, pppc[1])
    #elif pppc[0] == 'q':
    #    return newFormDedekindSumNaiveFastPrimePower(dChar1, dChar2, gamma, pppc[1][0], pppc[1][1])
    #TODO: CHANGE IF BELOW TO ELIF ABOVE ONCE PRIME CAPABILITES ARE RESTORED
    if pppc[0] == 'q':
        return newFormDedekindSumNaiveFastPrimePower(dChar1, dChar2, gamma, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return newFormDedekindSumNaiveFastComposite(dChar1, dChar2, gamma, q1 * q2)

#def newFormDedekindSumNaiveFastPrime(dChar1, dChar2, gamma, p):

def newFormDedekindSumNaiveFastPrimePower(dChar1, dChar2, gamma, p, k):
    n = p ** k
    reps = cosetRepsSLTwoZOverGammaZeroPrimePower(p, k)
    inHChecker = inGroupChecker(Gamma0(n))
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
    sum = dSp[letter][0]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
        sum = dSp[letter][0] + dSp[letter][1] * sum
    return sum

def newFormDedekindSumNaiveFastComposite(dChar1, dChar2, gamma, n):
    reps = cosetRepsSLTwoZOverGammaZeroComposite(n)
    inHChecker = inGroupChecker(Gamma0(n))
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
    sum = dSp[letter][0]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
        sum = dSp[letter][0] + dSp[letter][1] * sum
    return sum

def newFormDedekindSumNaiveFastGeneralized(dChar1, dChar2, gamma):
    n = modulus(dChar1) * modulus(dChar2)
    reps = cosetRepsSLTwoZOverGammaZeroGeneralized(n)
    inHChecker = inGroupChecker(Gamma0(n))
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
    sum = dSp[letter][0]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = [newFormDedekindSum(dChar1, dChar2, letter), psiChar(dChar1, dChar2, letter)]
        sum = dSp[letter][0] + dSp[letter][1] * sum
    return sum

def newFormDedekindSumNaiveFastGeneralizedPrecompute(dChar1, dChar2, gamma, chprPath):
    n = modulus(dChar1) * modulus(dChar2)
    reps = cosetRepsSLTwoZOverGammaZeroGeneralized(n)
    inHChecker = inGroupChecker(Gamma0(n))
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(gamma)), n)
    dSp = readAllChpr(chprPath)
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    msl = matrixString(letter)
    if not msl in dSp:
        print(msl)
        dSp[msl] = newFormDedekindSum(dChar1, dChar2, letter)
    sum = dSp[msl]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        msl = matrixString(letter)
        if not msl in dSp:
            print(msl)
            dSp[msl] = newFormDedekindSum(dChar1, dChar2, letter)
        sum = dSp[msl] + psiChar(dChar1, dChar2, letter) * sum
    return sum

def newFormDedekindSumFast(dChar1, dChar2, gamma):
    # computes the new form dedekind sum of gamma given two primative characters with similar parity
    # uses fast technique developed during 22 REU
    q1, q2 = modulus(dChar1), modulus(dChar2)
    pppc = pOppOc(q1 * q2)
    #if pppc[0] == 'p':
    #    return newFormDedekindSumFastPrime(dChar1, dChar2, gamma, pppc[1])
    #elif pppc[0] == 'q':
    #    return newFormDedekindSumFastPrimePower(dChar1, dChar2, gamma, pppc[1][0], pppc[1][1])
    #TODO: CHANGE IF BELOW TO ELIF ABOVE ONCE PRIME CAPABILITES ARE RESTORED
    if pppc[0] == 'q':
        return newFormDedekindSumFastPrimePower(dChar1, dChar2, gamma, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return newFormDedekindSumFastComposite(dChar1, dChar2, gamma, q1 * q2)

#def newFormDedekindSumFastPrime(dChar1, dChar2, gamma, p):

def newFormDedekindSumFastPrimePower(dChar1, dChar2, gamma, p, k):
    n = p ** k
    G0G1reps = cosetRepsGammaZeroOverGammaOnePrimePower(p, k)
    reps = cosetRepsSLTwoZOverGammaOnePrimePower(p, k)
    inHChecker = inGroupChecker(Gamma1(n))
    G0G1 = findCosetRepsLeft(G0G1reps, inHChecker, gamma, n)
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(G0G1 ** (-1) * gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
    sum = dSp[letter]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
        sum += dSp[letter]
    return newFormDedekindSum(dChar1, dChar2, G0G1) + psiChar(dChar1, dChar2, G0G1) * sum

def newFormDedekindSumFastComposite(dChar1, dChar2, gamma, n):
    G0G1reps = cosetRepsGammaZeroOverGammaOneComposite(n)
    reps = cosetRepsSLTwoZOverGammaOneComposite(n)
    inHChecker = inGroupChecker(Gamma1(n))
    G0G1 = findCosetRepsLeft(G0G1reps, inHChecker, gamma, n)
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(G0G1 ** (-1) * gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
    sum = dSp[letter]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
        sum += dSp[letter]
    return newFormDedekindSum(dChar1, dChar2, G0G1) + psiChar(dChar1, dChar2, G0G1) * sum

def newFormDedekindSumFastGeneralized(dChar1, dChar2, gamma):
    n = modulus(dChar1) * modulus(dChar2)
    G0G1reps = cosetRepsGammaZeroOverGammaOneGeneralized(n)
    reps = cosetRepsSLTwoZOverGammaOneGeneralized(n)
    inHChecker = inGroupChecker(Gamma1(n))
    G0G1 = findCosetRepsLeft(G0G1reps, inHChecker, gamma, n)
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(G0G1 ** (-1) * gamma)), n)
    dSp = {}
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    if not letter in dSp:
        dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
    sum = dSp[letter]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        if not letter in dSp:
            dSp[letter] = newFormDedekindSum(dChar1, dChar2, letter)
        sum += dSp[letter]
    return newFormDedekindSum(dChar1, dChar2, G0G1) + psiChar(dChar1, dChar2, G0G1) * sum

def newFormDedekindSumFastGeneralizedPrecompute(dChar1, dChar2, gamma, chprPath):
    n = modulus(dChar1) * modulus(dChar2)
    G0G1reps = cosetRepsGammaZeroOverGammaOneGeneralized(n)
    reps = cosetRepsSLTwoZOverGammaOneGeneralized(n)
    inHChecker = inGroupChecker(Gamma1(n))
    G0G1 = findCosetRepsLeft(G0G1reps, inHChecker, gamma, n)
    word = reidemeisterRewriteReps(reps, inHChecker, TSDecompToRewritingTape(TSDecomp(G0G1 ** (-1) * gamma)), n)
    dSp = readAllChpr(chprPath)
    l = len(word)
    letter = word[-1]
    letter.set_immutable()
    msl = matrixString(letter)
    if not msl in dSp:
        print(msl)
        dSp[msl] = newFormDedekindSum(dChar1, dChar2, letter)
    sum = dSp[msl]
    for i in tqdm(range(l - 1)):
        letter = word[l - 2 - i]
        letter.set_immutable()
        msl = matrixString(letter)
        if not msl in dSp:
            nfds = newFormDedekindSum(dChar1, dChar2, letter)
            print(msl, nfds)
            dSp[msl] = nfds
        sum += dSp[msl]
    return newFormDedekindSum(dChar1, dChar2, G0G1) + psiChar(dChar1, dChar2, G0G1) * sum

#def cdMaxNorm(norm = 10000):
#    pairs = []
#    for c in range(ceil(-1 * sqrt(norm)), floor(sqrt(norm))):
#        for d in range(ceil(-1 * sqrt(norm - c ** 2)), floor(sqrt(norm - c ** 2))):
#            if gcd(c, d) == 1:
#                pairs.append([c, d])
#    return pairs

#def newformEisensteinSeries(dChar1, dChar2, z, s, norm = 10000):
#    #TODO: needs fixing
#    if isEven(dChar1) == isEven(dChar2):
#        q1, q2 = modulus(dChar1), modulus(dChar2)
#        sum = 0
#        for pair in cdMaxNorm(norm):
#            c, d = pair[0], pair[1]
#            sum += ((q2 * z.imag()) ** s * dChar1(c) * dChar2(d)) / (c * q2 * z + d) ** (2 * s)
#        return sum / 2
#    return None

#def completedEisensteinSeries(dChar1, dChar2, z, s, norm = 10000):
#    #TODO: needs fixing
#    if isEven(dChar1) == isEven(dChar2):
#        q1, q2 = modulus(dChar1), modulus(dChar2)
#        return ((q2 / pi) ** s / gaussSum(dChar2)) * gammaFunction(s) * LFunction(2 * s, prodDCharacters(dChar1, dChar2)) * newformEisensteinSeries(dChar1, dChar2, z, s, norm)
#    return None