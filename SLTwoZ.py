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

def TSDecompToRewritingTape(TS):
    #takes a TS decomposition and creates a rewriting tape for the reidemeister schreier rewrititng process
    tape = [[matrix.identity(2), S ** (2 * TS[-1])]]
    for i in range(len(TS) - 2):
        tape.append([tape[-1][0] * tape[-1][1], T ** TS[i]])
        tape.append([tape[-1][0] * tape[-1][1], S])
    tape.append([tape[-1][0] * tape[-1][1], T ** TS[-2]])
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

def schreierLemmaReps(G, reps, inHChecker, n, verbose = True):
    #given G and a set of coset representatives of G/H and a method to check if an element is in H, this computes the generators of H
    gType = congSubGroupType(G)
    if gType == '0':
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

#def reidemeisterRewrite():

#def isSchreierTransversal(transversal):

def reidemeisterSchreierRewrite(G, H, rewritingTape, p):
    #computes a rewritten representation of a word in G as a word in H
    return [U(G, H, findCoset(G, H, stage[0], p), stage[1], p) for stage in rewritingTape]

def reidemeisterSchreierRewriteReps(reps, inHChecker, rewritingTape, p):
    #computes a rewritten representation of a word in G as a word in H given coset representatives G/H and a method to check if an element is in H
    return [UReps(reps, inHChecker, findCosetReps(reps, inHChecker, stage[0], p), stage[1], p) for stage in rewritingTape]

###############################
###############################
###                         ###
###   COSET INVESTIGATION   ###
###                         ###
###############################
###############################

def cosetRepsSLTwoZOverGammaZeroPrime(p):
    #returns a set of coset representatives of SL2Z/Gamma0(p), these representatives form a Schreier transversal
    reps = [matrix.identity(2)]
    for i in range(0, p):
        reps.append(S * T ** i)
    return reps

def cosetRepsSLTwoZOverGammaZeroPrimePower(p, k):
    #returns a set of coset representatives of SL2Z/Gamma0(p^k), these representatives for a Schreier transversal
    reps = cosetRepsSLTwoZOverGammaZeroPrime(p)
    for i in range(1, p ** (k - 1)):
        reps.append(S * T ** (i * p) ** S)
    return reps

#def cosetRepsSLTwoZOverGammaZeroComposite(n):
#    #returns a set of coset representatives of SL2Z/Gamma0(n), these representatives for a Schreier transversal (?)

#def cosetRepsSLTwoZOverGammaZero(n):
#    #returns a set of coset representatives of SL2Z/Gamma0(n)

def cosetRepsSLTwoZOverGammaOnePrime(p):
    #returns a set of coset representatives of SL2Z/Gamma1(p)
    return doubleQuotientLemma(cosetRepsSLTwoZOverGammaZeroPrime(p), cosetRepsGammaZeroOverGammaOnePrime(p))

def cosetRepsSLTwoZOverGammaOnePrimePower(p, k):
    #returns a set of coset representatives of SL2Z/Gamma1(p^k)
    return doubleQuotientLemma(cosetRepsSLTwoZOverGammaZeroPrimePower(p, k), cosetRepsGammaZeroOverGammaOnePrimePower(p, k))

#def cosetRepsSLTwoZOverGammaOneComposite(n):
#    #returns a set of coset representatives of SL2Z/Gamma1(n)

#def cosetRepsSLTwoZOverGammaOne(n):
#    #returns a set of coset representatives of SL2Z/Gamma1(n)

def cosetRepsGammaZeroOverGammaOnePrime(p):
    #returns a set of coset representatives of Gamma0(p)/Gamma1(p)
    reps = []
    for i in range(1, p):
        if gcd(i, p) == 1:
            ii = inverse_mod(i, p)
            reps.append(matrix(ZZ, [[i, (i * ii - 1) // p], [p, ii]]))
    return reps

def cosetRepsGammaZeroOverGammaOnePrimePower(p, k):
    #returns a set of coset representatives of Gamma0(p^k)/Gamma1(p^k)
    return cosetRepsGammaZeroOverGammaOnePrime(p ** k)

#def cosetRepsGammaZeroOverGammaOneComposite(n):
#    #returns a set of coset representatives of Gamma0(n)/Gamma1(n)

#def cosetRepsGammaZeroOverGammaOne(n):
#    #returns a set of coset representatives of Gamma0(n)/Gamma1(n)

def cosetRepsPrime(G, H, p):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZeroPrime(p)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOnePrime(p)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOnePrime(p)

def cosetRepsPrimePower(G, H, p, k):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(p^k), or Gamma1(p^k)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZeroPrimePower(p, k)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOnePrimePower(p, k)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOnePrimePower(p, k)

#def cosetRepsComposite(G, H, n):
#    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(n), or Gamma1(n)

#def cosetReps(G, H, n):
#    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(n), or Gamma1(n)

def findCosetSLTwoZOverGammaZeroPrime(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(p)
    reps = cosetRepsSLTwoZOverGammaZeroPrime(p)
    inHChecker = inGroupChecker(Gamma0(p))
    return findCosetReps(reps, inHChecker, element, p)

def findCosetSLTwoZOverGammaZeroPrimePower(element, p, k):
    # given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(p^k)
    reps = cosetRepsSLTwoZOverGammaZeroPrimePower(p, k)
    inHChecker = inGroupChecker(Gamma0(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

#def findCosetSLTwoZOverGammaZeroComposite(element, n):
#    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(n)

#def findCosetSLTwoZOerGammaZero(element, n):
#    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(n)

def findCosetSLTwoZOverGammaOnePrime(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(p)
    reps = cosetRepsSLTwoZOverGammaOnePrime(p)
    inHChecker = inGroupChecker(Gamma1(p))
    return findCosetReps(reps, inHChecker, element, p)

def findCosetSLTwoZOverGammaOnePrimePower(element, p, k):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(p^k)
    reps = cosetRepsSLTwoZOverGammaOnePrimePower(p, k)
    inHChecker = inGroupChecker(Gamma1(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

#def findCosetSLTwoZOverGammaOneComposite(element, n):
#    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(n)


#def findCosetSLTwoZOverGammaOne(element, n):
#    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(n)


def findCosetGammaZeroOverGammaOnePrime(element, p):
    #given an element in Gamma0(p), this computes the coset representative of the element in Gamma0(p)/Gamma1(p)
    reps = cosetRepsGammaZeroOverGammaOnePrime(p)
    inHChecker = inGroupChecker(Gamma1(p))
    return findCosetReps(reps, inHChecker, element, p)

def findCosetGammaZeroOverGammaOnePrimePower(element, p, k):
    #given an element in Gamma0(p^k), this computes the coset representative of the element in Gamma0(p^k)/Gamma1(p^k)
    reps = cosetRepsGammaZeroOverGammaOnePrimePower(p, k)
    inHChecker = inGroupChecker(Gamma1(p ** k))
    return findCosetReps(reps, inHChecker, element, p ** k)

#def findCosetGammaZeroOverGammaOneComposite(element, n):
#    #given an element in Gamma0(n), this computes the coset representative of the element in Gamma0(n)/Gamma1(n)


#def findCosetGammaZeroOverGammaOne(element, n):
#    #given an element in Gamma0(n), this computes the coset representative of the element in Gamma0(n)/Gamma1(n)

def findCosetPrime(G, H, element, p):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetSLTwoZOverGammaZeroPrime(element, p)
    elif gType == 'S' and hType == '1':
        return findCosetSLTwoZOverGammaOnePrime(element, p)
    elif gType == '0' and hType == '1':
        return findCosetGammaZeroOverGammaOnePrime(element, p)

def findCosetPrimePower(G, H, element, p, k):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(p^k), or Gamma1(p^k)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetsSLTwoZOverGammaZeroPrimePower(element, p, k)
    if gType == 'S' and hType == '1':
        return findCosetsSLTWoZOverGammaOnePrimePower(element, p, k)
    if gType == '0' and hType == '1':
        return findCosetsGammaZeroOverGammaOnePrimePower(element, p, k)

#def findCosetComposite(G, H, element, n):

#def findCoset(G, H, element, n):

def findCosetReps(reps, inHChecker, element, n):
    #given a set of coset representatives of G/H and a method to check if an element is in H, this computes the coset representative of the element in G/H
    for rep in reps:
        if inHChecker(element * rep ** (-1), n):
            return rep

######################################
######################################
###                                ###
###   GROUP ACTION INVESTIGATION   ###
###                                ###
######################################
######################################

#def SLTwoZOverGammaZeroGroupAction(p, groupAction = 'S'):
#    #describes the group action of the generators of SL2Z on the cosets of SL2Z/Gamma0(p)
#    if groupAction == 'S':
#        outDict = {-1: 0, 0: -1}
#        for i in range(1, p):
#            outDict[i] = (-1 * pow(i, p - 2, p)) % p
#        return outDict
#    elif groupAction == 'T':
#        outDict = {-1: -1}
#        for i in range(p):
#            outDict[i] = (i + 1) % p
#        return outDict
#    else:
#        print("invalid group action: not a generator of SL(2,ZZ)")
#        return -1

#def SLTwoZOverGammaOneGroupAction(p, groupAction = ''):

#def GammaZeroOverGammaOneGroupAction(p, groupAction = ''):