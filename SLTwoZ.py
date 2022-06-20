from base import *

#R = matrix(ZZ, [[0, -1], [1, 1]])
S = matrix(ZZ, [[0, -1], [1, 0]])
T = matrix(ZZ, [[1, 1], [0, 1]])
#U = matrix(ZZ, [[1, 0], [1, 1]])

#def V(p):
#    pr = primitiveRoot(p)
#    a, b = pr, pow(pr, p - 2, p)
#    return S * T ** a * S * T ** b * S * T ** a, a

#def R(mu, v, p):
#    muStar = (mu * v ** 2 - v) % p
#    vStar = (-1 * pow(v, p - 2, p)) % p
#    return T ** mu * S * T ** v * S * T ** (-1 * vStar) * S * T ** (-1 * muStar)

#def gammaPrimeGen(l, mu, v, p):
#    VMatrix, a = V(p)
#    s = Mod(v, p).log(Mod(a, p))
#    lStar = (l + s) % ((p - 1) // 2)
#    return VMatrix ** l * R(mu, v, p) * VMatrix ** (-lStar)

#def gammaPrimeGens(p):
#    return [gammaPrimeGen(l, mu, v, p) for v in range(1, p) for mu in range(p) for l in range((p - 3) // 2 + 1)]

def rotateMatrixCW(m):
    return m.transpose() * matrix(ZZ, [[0, 1], [1, 0]])

def rotateMatrixCCW(m):
    return (m * matrix(ZZ, [[0, 1], [1, 0]])).transpose()

def contFrac(a, c):
    out = []
    while c != 0:
        out.append(ceil(a / c))
        a, c = c, (-1 * a) % c
    return out

def contFracTwos(p):
  Ns = {}
  for k in range(p - 1):
    for n in range(2, p + 1):
      if (n * (k + 1) - k) % p == 0:
        Ns[k] = n
  return Ns

def buildMatrix(a, c):
    m = matrix.identity(2)
    for pow in contFrac(a, c):
        m = m * T ** pow * S
    return m

def computeOrder(m):
    for k in range(2, 7):
        if m ** k == matrix.identity(2):
            return k
    return -1

def inSLTwoZ(m):
    return True if m.determinant() == 1 else False

def inGammaZero(m, n):
    return inSLTwoZ(m) and (m[1][0] % n == 0)

def gammaZeroPrimeGens(p):
    gens = [-1 * matrix.identity(2), T ** 1, T ** (-1)]
    cFTs = contFracTwos(p)
    for k in range(p - 1):
        gens.append(S * T ** cFTs[k] * S * (T ** 2 * S) ** k)
    return gens

def inGammaOne(m, n):
    return inGammaZero(m, n) and (m[0][0] % n == 1) and (m[1][1] % n == 1)

def inGamma(m, n):
    return inGammaZero(m, n) and (m[0][1] % n == 1)

def SLTwoZOverGammaZeroGroupAction(p, groupAction = 's'):
    if groupAction == 's' or groupAction == 'S':
        outDict = {-1: 0, 0: -1}
        for i in range(1, p):
            outDict[i] = (-1 * pow(i, p - 2, p)) % p
        return outDict
    elif groupAction == 't' or groupAction == 'T':
        outDict = {-1: -1}
        for i in range(p):
            outDict[i] = (i + 1) % p
        return outDict
    else:
        print("invalid group action: not a generator of SL(2,ZZ)")
        return -1

#ef SLTwoZOverGammaOneGroupAction(p, groupAction = ''):

#def GammaZeroOverGammaOneGroupAction(p, groupAction = ''):