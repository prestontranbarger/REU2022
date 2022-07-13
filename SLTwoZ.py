from base import *

S = matrix(ZZ, [[0, -1], [1, 0]])
T = matrix(ZZ, [[1, 1], [0, 1]])

S.set_immutable()
T.set_immutable()

def matrixMod(m, p):
    #returns a matrix modulo p element wise
    return matrix(ZZ, [[m[0][0] % p, m[0][1] % p], [m[1][0] % p, m[1][1] % p]])

def rotateMatrixCW(m):
    #rotates a matrix's entries clockwise
    return m.transpose() * matrix(ZZ, [[0, 1], [1, 0]])

def rotateMatrixCCW(m):
    #rotates a matrix's entries counter-clockwise
    return (m * matrix(ZZ, [[0, 1], [1, 0]])).transpose()

def computeOrder(m):
    #computes the order of a matrix in SL2Z
    for k in range(2, 7):
        if m ** k == matrix.identity(2):
            return k
    return -1

def contFrac(a, c):
    #returns the continuted fraction expansion of a/c
    if c < 0:
        a, c = -1 * a, -1 * c
    out = []
    while c != 0:
        out.append(ceil(a / c))
        a, c = c, (-1 * a) % c
    return out

def TSDecomp(m):
    #returns the TS decomposition of a matrix
    TS = []
    while m[1][0] != 0:
        exp = -1 * floor(m[0][0] / m[1][0])
        m = S * T ** exp * m
        TS.append(-1 * exp)
    pm = int(m[0][0] == 1)
    TS.append((2 * pm - 1) * m[0][1])
    TS.append((len(TS) + pm) % 2)
    return TS

def TSDecompToRewritingTape(tsd):
    #takes a TS decomposition and creates a rewriting tape for the reidemeister schreier rewrititng process
    tape = [[matrix.identity(2), S ** (2 * tsd[-1])]]
    for i in range(len(tsd) - 2):
        tape.append([tape[-1][0] * tape[-1][1], T ** tsd[i]])
        tape.append([tape[-1][0] * tape[-1][1], S])
    tape.append([tape[-1][0] * tape[-1][1], T ** tsd[-2]])
    return tape

def buildMatrix(a, c):
    #builds a matrix with left column [a, c]
    m = matrix.identity(2)
    for pow in contFrac(a, c):
        m = m * T ** pow * S
    return m

def buildMatrixFromCF(cF):
    #builds a matrix from a specific continued fraction expansion of a/c
    m = matrix.identity(2)
    for pow in cF:
        m = m * T ** pow * S
    return m

def buildMatrixFromTS(TS):
    #builds a matrix from a TS decomposition
    m = S ** (2 * TS[-1])
    for pow in TS[:-1]:
        m = m * T ** pow * S
    return m * S ** (-1)

def congSubGroupType(csg):
    #checks what congruence subgroup an object belongs to
    #works for all p, note that Gamma0(2) and Gamma1(2) are the same
    if isinstance(csg, type(Gamma1(3))):
        return '1'
    elif isinstance(csg, type(Gamma0(3))):
        return '0'
    return 'S'

def inSLTwoZ(m):
    #checks if m is in sl2z
    return True if m.determinant() == 1 else False

def inGammaZero(m, n):
    #checks if m is in Gamma0(n)
    return inSLTwoZ(m) and (m[1][0] % n == 0)

def inGammaOne(m, n):
    #checks if m is in Gamma1(n)
    return inGammaZero(m, n) and (m[0][0] % n == 1) and (m[1][1] % n == 1)

def inGamma(m, n):
    #checks if m is in Gamma(n)
    return inGammaOne(m, n) and (m[0][1] % n == 0)

def inGroupChecker(G):
    #returns the method which checks if an element belongs to the group corresponding to the parameter
    if congSubGroupType(G) == '1':
        return inGammaOne
    elif congSubGroupType(G) == '0':
        return inGammaZero
    return inSLTwoZ

#def gammaZeroPrimeGens(p):
#    gens = [-1 * matrix.identity(2), T ** 1, T ** (-1)]
#    cFTs = contFracTwos(p)
#    for k in range(p - 1):
#        gens.append(S * T ** cFTs[k] * S * (T ** 2 * S) ** k)
#    return gens

############################
############################
###                      ###
###   IMPORTANT LEMMAS   ###
###                      ###
############################
############################

def doubleQuotientLemma(reps1, reps2):
    #given coset representatives of G/H and H/K, computes coset representatives of G/K
    return [rep1 * rep2 for rep1 in reps1 for rep2 in reps2]

def U(G, H, t, s, n):
    #given quotient group G/H, this computes the schreier transformation of t, s
    return t * s * findCoset(G, H, t * s, n) ** -1

def UReps(reps, inHChecker, t, s, n):
    #given a set of coset representatives of G/H and a method to check if an element is in H, this computes the schreier transformation of t, s
    return t * s * findCosetReps(reps, inHChecker, t * s, n) ** -1

def schreierLemma(G, H, n, verbose = False):
    #given quotient group G/H, this computes the generators of H
    reps = cosetReps(G, H, n)
    gType = congSubGroupType(G)
    if gType == '0':
        gensIn = schreierLemma(SL2Z, Gamma0(n), n)
    elif gType == 'S':
        gensIn = [S, T]
    if verbose:
        gensOut = {}
        for rep in reps:
            for genIn in gensIn:
                genOut = U(G, H, rep, genIn, n)
                rep.set_immutable()
                genOut.set_immutable()
                gensOut[(rep, genIn)] = genOut
        return gensOut
    else:
        gensOut = []
        for rep in reps:
            for genIn in gensIn:
                genOut = U(G, H, rep, genIn, n)
                genOut.set_immutable()
                gensOut.append(genOut)
        return list(set(gensOut))

def schreierLemmaReps(G, reps, inHChecker, n, verbose = False):
    #given G and a set of coset representatives of G/H and a method to check if an element is in H, this computes the generators of H
    gType = congSubGroupType(G)
    if gType == '0':
        #TODO: IMPROVE BY USING SCHREIER LEMMA REPS
        gensIn = schreierLemma(SL2Z, Gamma0(n), n)
    elif gType == 'S':
        gensIn = [S, T]
    if verbose:
        gensOut = {}
        for rep in reps:
            for genIn in gensIn:
                genOut = UReps(reps, inHChecker, rep, genIn, n)
                rep.set_immutable()
                genOut.set_immutable()
                gensOut[(rep, genIn)] = genOut
        return gensOut
    else:
        gensOut = []
        for rep in reps:
            for genIn in gensIn:
                genOut = UReps(reps, inHChecker, rep, genIn, n)
                genOut.set_immutable()
                gensOut.append(genOut)
        return list(set(gensOut))

#def isSchreierTransversal(transversal):

def reidemeisterRewrite(G, H, rewritingTape, n):
    #computes a rewritten representation of a word in G as a word in H
    return [U(G, H, findCoset(G, H, stage[0], n), stage[1], n) for stage in rewritingTape]

def reidemeisterRewriteReps(reps, inHChecker, rewritingTape, n):
    #computes a rewritten representation of a word in G as a word in H given coset representatives G/H and a method to check if an element is in H
    return [UReps(reps, inHChecker, findCosetReps(reps, inHChecker, stage[0], n), stage[1], n) for stage in rewritingTape]

def projectiveLift(c, d, N):
    if c != 0 and d != 0:
        while gcd(c, d) != 1:
            d += N
        b = (-1 * inverse_mod(c, d)) % d
        a = (1 + b * c) // d
        return matrix(ZZ, [[a, b],
                           [c, d]])
    elif c == 0:
        di = inverse_mod(d, N)
        return matrix(ZZ, [[di, (d * di - 1) // N],
                           [N, d]])
    elif d == 0:
        ci = (-1 * inverse_mod(c, N)) % N
        return matrix(ZZ, [[(c * ci + 1) // N, ci],
                           [c, N]])

def PRing(n):
    elements = []
    for c in range(n):
        for d in range(n):
            if gcd(gcd(c, d), n) == 1:
                elements.append((c, d))
    return elements

def POneRing(n):
    units = []
    for i in range(n):
        if gcd(i, n) == 1:
            units.append(i)
    pairs = PRing(n)
    for pair in pairs:
        for unit in units:
            new = (pair[0] * unit % n, pair[1] * unit % n)
            if new in pairs and unit != 1:
                pairs.remove(new)
    return pairs

###############################
###############################
###                         ###
###   COSET INVESTIGATION   ###
###                         ###
###############################
###############################

#cosets SL2Z/Gamma0(n)

def cosetRepsSLTwoZOverGammaZeroPrime(p):
    #returns a set of coset representatives of SL2Z/Gamma0(p)
    reps = [matrix.identity(2)]
    for i in range(0, p):
        reps.append(S * T ** i)
    return reps

def findCosetSLTwoZOverGammaZeroPrime(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(p)
    reps = cosetRepsSLTwoZOverGammaZeroPrime(p)
    inHChecker = inGroupChecker(Gamma0(p))
    return findCosetReps(reps, inHChecker, element, p)

def cosetRepsSLTwoZOverGammaZeroPrimePower(p, k):
    #returns a set of coset representatives of SL2Z/Gamma0(p^k)
    reps = cosetRepsSLTwoZOverGammaZeroPrime(p ** k)
    for i in range(1, p ** (k - 1)):
        reps.append(S * T ** (i * p) * S)
    return reps

def findCosetSLTwoZOverGammaZeroPrimePower(element, p, k):
    # given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(p^k)
    reps = cosetRepsSLTwoZOverGammaZeroPrimePower(p, k)
    inHChecker = inGroupChecker(Gamma0(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

def cosetRepsSLTwoZOverGammaZeroComposite(n):
    #returns a set of coset representatives of SL2Z/Gamma0(n)
    return [projectiveLift(pair[0], pair[1], n) for pair in POneRing(n)]

def findCosetSLTwoZOverGammaZeroComposite(element, n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(n)
    reps = cosetRepsSLTwoZOverGammaZeroComposite(n)
    inHChecker = inGroupChecker(Gamma0(n))
    return findCosetReps(reps, inHChecker, element, n)

def cosetRepsSLTwoZOverGammaZero(n):
    #returns a set of coset representatives of SL2Z/Gamma0(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return cosetRepsSLTwoZOverGammaZeroPrime(pppc[1])
    elif pppc[0] == 'q':
        return cosetRepsSLTwoZOverGammaZeroPrimePower(pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return cosetRepsSLTwoZOverGammaZeroComposite(n)

def findCosetSLTwoZOverGammaZero(element, n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return findCosetSLTwoZOverGammaZeroPrime(element, pppc[1])
    elif pppc[0] == 'q':
        return findCosetSLTwoZOverGammaZeroPrimePower(element, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return findCosetSLTwoZOverGammaZeroComposite(element, n)

def cosetRepsSLTwoZOverGammaZeroGeneralized(n):
    #returns a set of coset representatives of SL2Z/Gamma0(n), uses the projective lift regardless of n
    #designed to reduce total method calls
    return [projectiveLift(pair[0], pair[1], n) for pair in POneRing(n)]

def findCosetSLTwoZOverGammaZeroGeneralized(element, n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(n), uses the projective lift regardless of n
    #designed to reduce total method calls
    reps = cosetRepsSLTwoZOverGammaZeroGeneralized(n)
    inHChecker = inGroupChecker(Gamma0(n))
    return findCosetReps(reps, inHChecker, element, n)

#cosets SL2Z/Gamma1(n)

def cosetRepsSLTwoZOverGammaOnePrime(p):
    #returns a set of coset representatives of SL2Z/Gamma1(p)
    return doubleQuotientLemma(cosetRepsSLTwoZOverGammaZeroPrime(p), cosetRepsGammaZeroOverGammaOnePrime(p))

def findCosetSLTwoZOverGammaOnePrime(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(p)
    reps = cosetRepsSLTwoZOverGammaOnePrime(p)
    inHChecker = inGroupChecker(Gamma1(p))
    return findCosetReps(reps, inHChecker, element, p)

def cosetRepsSLTwoZOverGammaOnePrimePower(p, k):
    #returns a set of coset representatives of SL2Z/Gamma1(p^k)
    return doubleQuotientLemma(cosetRepsSLTwoZOverGammaZeroPrimePower(p, k), cosetRepsGammaZeroOverGammaOnePrimePower(p, k))

def findCosetSLTwoZOverGammaOnePrimePower(element, p, k):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(p^k)
    reps = cosetRepsSLTwoZOverGammaOnePrimePower(p, k)
    inHChecker = inGroupChecker(Gamma1(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

def cosetRepsSLTwoZOverGammaOneComposite(n):
    #returns a set of coset representatives of SL2Z/Gamma1(n)
    return [projectiveLift(pair[0], pair[1], n) for pair in PRing(n)]

def findCosetSLTwoZOverGammaOneComposite(element, n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(n)
    reps = cosetRepsSLTwoZOverGammaOneComposite(n)
    inHChecker = inGroupChecker(Gamma1(n))
    return findCosetReps(reps, inHChecker, element, n)

def cosetRepsSLTwoZOverGammaOne(n):
    #returns a set of coset representatives of SL2Z/Gamma1(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return cosetRepsSLTwoZOverGammaOnePrime(pppc[1])
    elif pppc[0] == 'q':
        return cosetRepsSLTwoZOverGammaOnePrimePower(pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return cosetRepsSLTwoZOverGammaOneComposite(n)

def findCosetSLTwoZOverGammaOne(element, n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return findCosetSLTwoZOverGammaOnePrime(element, pppc[1])
    elif pppc[0] == 'q':
        return findCosetSLTwoZOverGammaOnePrimePower(element, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return findCosetSLTwoZOverGammaOneComposite(element, n)

def cosetRepsSLTwoZOverGammaOneGeneralized(n):
    #returns a set of coset representatives of SL2Z/Gamma1(n), uses the projective lift for all n
    #designed to reduce total method calls
    return [projectiveLift(pair[0], pair[1], n) for pair in PRing(n)]

def findCosetRepsSLTwoZOverGammaOneGeneralized(n):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(n), uses the projective lift for all n
    #designed to reduce total method calls
    reps = cosetRepsSLTwoZOverGammaOneGeneralized(n)
    inHChecker = inGroupChecker(Gamma1(n))
    return findCosetReps(reps, inHChecker, element, n)

#cosets Gamma0(n)/Gamma1(n)

def cosetRepsGammaZeroOverGammaOnePrime(p):
    #returns a set of coset representatives of Gamma0(p)/Gamma1(p)
    reps = []
    for i in range(1, p):
        if gcd(i, p) == 1:
            ii = inverse_mod(i, p)
            reps.append(matrix(ZZ, [[i, (i * ii - 1) // p], [p, ii]]))
    return reps

def findCosetGammaZeroOverGammaOnePrime(element, p):
    #given an element in Gamma0(p), this computes the coset representative of the element in Gamma0(p)/Gamma1(p)
    reps = cosetRepsGammaZeroOverGammaOnePrime(p)
    inHChecker = inGroupChecker(Gamma1(p))
    return findCosetReps(reps, inHChecker, element, p)

def cosetRepsGammaZeroOverGammaOnePrimePower(p, k):
    #returns a set of coset representatives of Gamma0(p^k)/Gamma1(p^k)
    return cosetRepsGammaZeroOverGammaOnePrime(p ** k)

def findCosetGammaZeroOverGammaOnePrimePower(element, p, k):
    #given an element in Gamma0(p^k), this computes the coset representative of the element in Gamma0(p^k)/Gamma1(p^k)
    reps = cosetRepsGammaZeroOverGammaOnePrimePower(p, k)
    inHChecker = inGroupChecker(Gamma1(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

def cosetRepsGammaZeroOverGammaOneComposite(n):
    #returns a set of coset representatives of Gamma0(n)/Gamma1(n)
    return cosetRepsGammaZeroOverGammaOnePrime(n)

def findCosetGammaZeroOverGammaOneComposite(element, n):
    #given an element in Gamma0(n), this computes the coset representative of the element in Gamma0(n)/Gamma1(n)
    reps = cosetRepsGammaZeroOverGammaOneComposite(n)
    inHChecker = inGroupChecker(Gamma1(n))
    return findCosetReps(reps, inHChecker, element, n)

def cosetRepsGammaZeroOverGammaOne(n):
    #returns a set of coset representatives of Gamma0(n)/Gamma1(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return cosetRepsGammaZeroOverGammaOnePrime(pppc[1])
    elif pppc[0] == 'q':
        return cosetRepsGammaZeroOverGammaOnePrimePower(pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return cosetRepsGammaZeroOverGammaOneComposite(n)

def findCosetGammaZeroOverGammaOne(element, n):
    #given an element in Gamma0(n), this computes the coset representative of the element in Gamma0(n)/Gamma1(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return findCosetGammaZeroOverGammaOnePrime(element, pppc[1])
    elif pppc[0] == 'q':
        return findCosetGammaZeroOverGammaOnePrimePower(element, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return findCosetGammaZeroOverGammaOneComposite(element, n)

def cosetRepsGammaZeroOverGammaOneGeneralized(n):
    #returns a set of coset representatives of Gamma0(n)/Gamma1(n)
    #designed to reduce total method calls
    reps = []
    for i in range(1, n):
        if gcd(i, n) == 1:
            ii = inverse_mod(i, n)
            reps.append(matrix(ZZ, [[i, (i * ii - 1) // n], [n, ii]]))
    return reps

def findCosetGammaZeroOverGammaOneGeneralized(element, n):
    #given an element in Gamma0(n), this computes the coset representative of the element in Gamma0(n)/Gamma1(n)
    #designed to reduce total method calls
    reps = cosetRepsGammaZeroOverGammaOneGeneralized(n)
    inHChecker = inGroupChecker(Gamma1(n))
    return findCosetReps(reps, inHChecker, element, n)

#general cosets

def cosetRepsPrime(G, H, p):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZeroPrime(p)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOnePrime(p)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOnePrime(p)

def findCosetPrime(G, H, element, p):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetSLTwoZOverGammaZeroPrime(element, p)
    elif gType == 'S' and hType == '1':
        return findCosetSLTwoZOverGammaOnePrime(element, p)
    elif gType == '0' and hType == '1':
        return findCosetGammaZeroOverGammaOnePrime(element, p)

def cosetRepsPrimePower(G, H, p, k):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(p^k), or Gamma1(p^k)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZeroPrimePower(p, k)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOnePrimePower(p, k)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOnePrimePower(p, k)

def findCosetPrimePower(G, H, element, p, k):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(p^k), or Gamma1(p^k)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetSLTwoZOverGammaZeroPrimePower(element, p, k)
    if gType == 'S' and hType == '1':
        return findCosetSLTwoZOverGammaOnePrimePower(element, p, k)
    if gType == '0' and hType == '1':
        return findCosetGammaZeroOverGammaOnePrimePower(element, p, k)

def cosetRepsComposite(G, H, n):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(n), or Gamma1(n)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZeroComposite(n)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOneComposite(n)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOneComposite(n)

def findCosetComposite(G, H, element, n):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(n), or Gamma1(n)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetSLTwoZOverGammaZeroComposite(element, p, k)
    elif gType == 'S' and hType == '1':
        return findCosetSLTwoZOverGammaOneComposite(element, n)
    elif gType == '0' and hType == '1':
        return findCosetGammaZeroOverGammaOneComposite(element, n)

def cosetReps(G, H, n):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(n), or Gamma1(n)
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return cosetRepsPrime(G, H, pppc[1])
    elif pppc[0] == 'q':
        return cosetRepsPrimePower(G, H, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return cosetRepsComposite(G, H, n)

def findCoset(G, H, element, n):
    pppc = pOppOc(n)
    if pppc[0] == 'p':
        return findCosetPrime(G, H, element, pppc[1])
    elif pppc[0] == 'q':
        return findCosetPrimePower(G, H, element, pppc[1][0], pppc[1][1])
    elif pppc[0] == 'n':
        return findCosetComposite(G, H, element, n)

def findCosetReps(reps, inHChecker, element, n):
    #given a set of coset representatives of G/H and a method to check if an element is in H, this computes the coset representative of the element in G/H
    for rep in reps:
        if inHChecker(element * rep ** (-1), n):
            return rep

def findCosetRepsLeft(reps, inHChecker, element, n):
    #given a set of LEFT coset representatives of G/H and a method to check if an element is in H, this computes the LEFT coset representative of the element in G/H
    for rep in reps:
        if inHChecker(rep ** (-1) * element, n):
            return rep

######################################
######################################
###                                ###
###   GROUP ACTION INVESTIGATION   ###
###                                ###
######################################
######################################

def SLTwoZOverGammaZeroGroupAction(p, groupAction = 'S'):
    #describes the group action of the generators of SL2Z on the cosets of SL2Z/Gamma0(p)
    if groupAction == 'S':
        outDict = {-1: 0, 0: -1}
        for i in range(1, p):
            outDict[i] = (-1 * pow(i, p - 2, p)) % p
        return outDict
    elif groupAction == 'T':
        outDict = {-1: -1}
        for i in range(p):
            outDict[i] = (i + 1) % p
        return outDict
    else:
        ("invalid group action: not a generator of SL(2,ZZ)")
        return -1

#def SLTwoZOverGammaOneGroupAction(p, groupAction = ''):

#def GammaZeroOverGammaOneGroupAction(p, groupAction = ''):