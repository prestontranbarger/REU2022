from base import *

R = matrix(ZZ, [[0, -1], [1, 1]])
S = matrix(ZZ, [[0, -1], [1, 0]])
T = matrix(ZZ, [[1, 1], [0, 1]])
U = matrix(ZZ, [[1, 0], [1, 1]])

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

def inGammaOne(m, n):
    return inGammaZero(m, n) and (m[0][0] % n == 1) and (m[1][1] % n == 1)

def inGamma(m, n):
    return inGammaZero(m, n) and (m[0][1] % n == 1)