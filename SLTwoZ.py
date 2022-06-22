from base import *

R = matrix(ZZ, [[0, -1], [1, 1]])
S = matrix(ZZ, [[0, -1], [1, 0]])
T = matrix(ZZ, [[1, 1], [0, 1]])
U = matrix(ZZ, [[1, 0], [1, 1]])

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
    TS.append(1 - int((-1) ** (len(TS) + pm - 2) == 1))
    return TS

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

def buildMatrixFromTS(tsDecomp):
    #builds a matrix from a TS decomposition
    m = S ** (2 * tsDecomp[-1])
    for pow in tsDecomp[:-1]:
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

#def gammaZeroPrimeGens(p):
#    gens = [-1 * matrix.identity(2), T ** 1, T ** (-1)]
#    cFTs = contFracTwos(p)
#    for k in range(p - 1):
#        gens.append(S * T ** cFTs[k] * S * (T ** 2 * S) ** k)
#    return gens

##########################
##                      ##
##   IMPORTANT LEMMAS   ##
##                      ##
##########################

def doubleQuotientLemma(reps1, reps2):
    #given coset representatives of G/H and H/K, computes coset representatives of G/K
    return [rep1 * rep2 for rep1 in reps1 for rep2 in reps2]

def schreierLemma(G, H, p):
    #given coset representatives of G/H and generators of G, computes generators of H
    reps = cosetReps(G, H, p)
    gType = congSubGroupType(G)
    if gType == '0':
        gensIn = schreierLemma(SL2Z, Gamma0(p), p)
    elif gType == 'S':
        gensIn = [S, T]
    gensOut = []
    for rep in reps:
        for genIn in gensIn:
            genOut = rep * genIn * findCoset(SL2Z, H, rep * genIn, p) ** -1
            genOut.set_immutable()
            gensOut.append(genOut)
    return list(set(gensOut))

#def reidemeisterRewrite():

#def isSchreierTransversal(transversal):

#def reidemeisterSchreierRewrite():

###########################
##                       ##
##   IMPORTANT METHODS   ##
##                       ##
###########################

def cosetRepsSLTwoZOverGammaZero(p):
    #returns a set of coset representatives of SL2Z/Gamma0(p), these representatives form a Schreier transversal
    reps = [matrix.identity(2)]
    for i in range(0, p):
        reps.append(S * T ** i)
    return reps

def cosetRepsGammaZeroOverGammaOne(p):
    #returns a set of coset representatives of SL2Z/Gamma1(p)
    reps = []
    for i in range(1, p):
        if gcd(i, p) == 1:
            ii = pow(i, p - 2, p)
            reps.append(matrix(ZZ, [[i, (i * ii - 1) // p], [p, ii]]))
    return reps

def cosetRepsSLTwoZOverGammaOne(p):
    #returns a set of coset representatives of Gamma0(p)/Gamma1(p)
    return doubleQuotientLemma(cosetRepsGammaZeroOverGammaOne(p), cosetRepsSLTwoZOverGammaZero(p))

def cosetReps(G, H, p):
    #returns a set of coset representatives of G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return cosetRepsSLTwoZOverGammaZero(p)
    elif gType == 'S' and hType == '1':
        return cosetRepsSLTwoZOverGammaOne(p)
    elif gType == '0' and hType == '1':
        return cosetRepsGammaZeroOverGammaOne(p)
    return [-1]

def findCosetSLTwoZOverGammaZero(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma0(p)
    reps = cosetRepsSLTwoZOverGammaZero(p)
    for rep in reps:
        if inGammaZero(element * rep ** (-1), p):
            return rep

def findCosetSLTwoZOverGammaOne(element, p):
    #given an element in SL2Z, this computes the coset representative of the element in SL2Z/Gamma1(p)
    reps = cosetRepsSLTwoZOverGammaOne(p)
    for rep in reps:
        if inGammaOne(element * rep ** (-1), p):
            return rep

def findCosetGammaZeroOverGammaOne(element, p):
    #given an element in Gamma0(p), this computes the coset representative of the element in Gamma0(p)/Gamma1(p)
    reps = cosetRepsGammaZeroOverGammaOne(p)
    for rep in reps:
        if inGammaOne(element * rep ** (-1), p):
            return rep

def findCoset(G, H, element, p):
    #given an element in G, this computes the coset representative of the element in G/H when H<G and G, H are any two of SL2Z, Gamma0(p), or Gamma1(p)
    gType, hType = congSubGroupType(G), congSubGroupType(H)
    if gType == 'S' and hType == '0':
        return findCosetSLTwoZOverGammaZero(element, p)
    elif gType == 'S' and hType == '1':
        return findCosetSLTwoZOverGammaOne(element, p)
    elif gType == '0' and hType == '1':
        return findCosetGammaZeroOverGammaOne(element, p)

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
        print("invalid group action: not a generator of SL(2,ZZ)")
        return -1

#def SLTwoZOverGammaOneGroupAction(p, groupAction = ''):

#def GammaZeroOverGammaOneGroupAction(p, groupAction = ''):