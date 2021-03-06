from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

import os
from os.path import isfile
from os.path import isdir

gensWritePath = "/precomputation"
charWritePath = "/precomputation"

def modPairs(N):
    pairs = []
    for j in range(3, N // 3 + 1):
        for k in range(3, N // j + 1):
            pairs.append((j, k))
    return pairs

def dCharString(dChar):
    l = str(dChar).split()
    out = l[3] + "c" + l[6] + ";"
    Tn = totient(modulus(dChar))
    zTn = CyclotomicField(Tn).gen()
    ns = [ss.split(" |--> ")[0].split()[-1] for ss in str(dChar).split(",")]
    for n in ns:
        p = 0
        for j in range(Tn):
            if zTn ** j == dChar(int(n)):
                p = j
        out += n + "-" + str(p) + "-"
    return out[:-1]

#def stringDChar(s):

# def precomputeGensGammaZero(n):
#     subPath = os.getcwd() + gensWritePath + str("/gammaZeroGens/")
#     filePath = subPath + "G0-" + str(n) + ".gens"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     if os.path.exists(filePath):
#         f = open(filePath, "r")
#         gens = [stringMatrix(line[:-1]) for line in f.readlines()]
#         f.close()
#     else:
#         reps = cosetRepsSLTwoZOverGammaZero(n)
#         inHChecker = inGroupChecker(Gamma0(n))
#         gens = schreierLemmaReps(SL2Z, reps, inHChecker, n)
#         f = open(filePath, "w")
#         f.writelines([matrixString(gen) + "\n" for gen in gens])
#         f.close()
#     return gens

# def precomputeGensGammaOne(n):
#     subPath = os.getcwd() + gensWritePath + str("/gammaOneGens/")
#     filePath = subPath + "G1-" + str(n) + ".gens"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     if os.path.exists(filePath):
#         f = open(filePath, "r")
#         gens = [stringMatrix(line[:-1]) for line in f.readlines()]
#         f.close()
#     else:
#         reps = cosetRepsSLTwoZOverGammaOne(n)
#         inHChecker = inGroupChecker(Gamma1(n))
#         gens = schreierLemmaReps(SL2Z, reps, inHChecker, n)
#         f = open(filePath, "w")
#         f.writelines([matrixString(gen) + "\n" for gen in gens])
#         f.close()
#     return gens

def precomputeMatrices(n):
    subPath = os.getcwd() + gensWritePath + str("/mtrxs/")
    filePath = subPath + str(n) + ".mtrxs"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        mtrxs = [stringMatrix(line[:-1]) for line in f.readlines()]
        f.close()
    else:
        repsDict, gRepsDict = cosetRepsSLTwoZOverGammaOneFast(n), cosetRepsGammaZeroOverGammaOneFast(n)
        mtrxs = []
        for key in gRepsDict:
            entry = gRepsDict[key]
            entry.set_immutable()
            mtrxs.append(entry)
        for key in repsDict:
            for exp in range(1, n + 1):
                tempMtrx = repsDict[key] * T ** exp
                entry = tempMtrx * repsDict[(tempMtrx[1][0] % n, tempMtrx[1][1] % n)] ** (-1)
                entry.set_immutable()
                mtrxs.append(entry)
            for exp in range(0, 3):
                tempMtrx = repsDict[key] * S ** exp
                entry = tempMtrx * repsDict[(tempMtrx[1][0] % n, tempMtrx[1][1] % n)] ** (-1)
                entry.set_immutable()
                mtrxs.append(entry)
        mtrxs = list(set(mtrxs))
        f = open(filePath, "w")
        f.writelines([matrixString(mtrx) + "\n" for mtrx in mtrxs])
        f.close()
    return mtrxs

# def precomputeCharacterPairs(dChar1, dChar2):
#     q1q2 = modulus(dChar1) * modulus(dChar2)
#     subPath, (G0Path, G1Path) = createCharacterPairFiles(dChar1, dChar2)
#     if G0Path != True:
#         f = open(G0Path, "w")
#         f.writelines([matrixString(gen) + complexString(newFormDedekindSum(dChar1, dChar2, gen)) + "\n" for gen in tqdm(precomputeGensGammaZero(q1q2))])
#         f.close()
#     if G1Path != True:
#         f = open(G1Path, "w")
#         f.writelines([matrixString(gen) + complexString(newFormDedekindSum(dChar1, dChar2, gen))+ "\n" for gen in tqdm(precomputeGensGammaOne(q1q2))])
#         f.close()

def precomputeCharacterPairs(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    subPath, filePath = createCharacterPairFile(dChar1, dChar2)
    if filePath != True:
        f = open(filePath, "w")
        f.writelines([matrixString(mtrx) + complexString(newFormDedekindSum(dChar1, dChar2, mtrx)) + "\n" for mtrx in tqdm(precomputeMatrices(n))])
        f.close()

# def createCharacterPairFiles(dChar1, dChar2):
#     q1q2 = modulus(dChar1) * modulus(dChar2)
#     str1, str2 = dCharString(dChar1), dCharString(dChar2)
#     subPath = os.getcwd() + gensWritePath + str("/characterPairs/")
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str1.split("c")[0] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str1.split(";")[0].split("c")[1] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str1.split(";")[1] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str2.split("c")[0] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str2.split(";")[0].split("c")[1] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     subPath += str2.split(";")[1] + "/"
#     if not isdir(subPath):
#         os.mkdir(subPath)
#     G0Path,\
#     G1Path = subPath + "G0-" + str(q1q2) + ".chpr",\
#              subPath + "G1-" + str(q1q2) + ".chpr"
#     if not os.path.exists(G0Path):
#         f = open(G0Path, "w")
#         f.close()
#     else:
#         G0Path = True
#     if not os.path.exists(G1Path):
#         f = open(G1Path, "w")
#         f.close()
#     else:
#         G1Path = True
#     return subPath, (G0Path, G1Path)

def createCharacterPairFile(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    str1, str2 = dCharString(dChar1), dCharString(dChar2)
    subPath = os.getcwd() + gensWritePath + str("/characterPairs/")
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split("c")[0] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split(";")[0].split("c")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str1.split(";")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split("c")[0] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split(";")[0].split("c")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    subPath += str2.split(";")[1] + "/"
    if not isdir(subPath):
        os.mkdir(subPath)
    filePath = subPath + str(n) + ".chpr"
    if not os.path.exists(filePath):
        f = open(filePath, 'w')
        f.close()
    else:
        filePath = True
    return subPath, filePath

# def chprPathFinder(dChar1, dChar2):
#     q1q2 = modulus(dChar1) * modulus(dChar2)
#     str1, str2 = dCharString(dChar1), dCharString(dChar2)
#     subPath = os.getcwd() + gensWritePath + \
#               str("/characterPairs/") + \
#               str1.split("c")[0] + "/" + \
#               str1.split(";")[0].split("c")[1] + "/" + \
#               str1.split(";")[1] + "/" + \
#               str2.split("c")[0] + "/" + \
#               str2.split(";")[0].split("c")[1] + "/" + \
#               str2.split(";")[1] + "/"
#     G0Path, \
#     G1Path = subPath + "G0-" + str(q1q2) + ".chpr", \
#              subPath + "G1-" + str(q1q2) + ".chpr"
#     return G0Path, G1Path

def chprPathFinder(dChar1, dChar2):
    n = modulus(dChar1) * modulus(dChar2)
    str1, str2 = dCharString(dChar1), dCharString(dChar2)
    subPath = os.getcwd() + gensWritePath + \
              str("/characterPairs/") + \
              str1.split("c")[0] + "/" + \
              str1.split(";")[0].split("c")[1] + "/" + \
              str1.split(";")[1] + "/" + \
              str2.split("c")[0] + "/" + \
              str2.split(";")[0].split("c")[1] + "/" + \
              str2.split(";")[1] + "/"
    filePath = subPath + str(n) + ".chpr"
    return filePath

# def readGens(path):
#     f = open(path, 'r')
#     gens = [line[:-1] for line in f.readlines()]
#     f.close()
#     return gens

def readMatrices(path):
    f = open(path, 'r')
    mtrxs = [line[:-1] for line in f.readlines()]
    f.close()
    return mtrxs

def readAllChpr(path):
    dict = {}
    f = open(path, 'r')
    for line in f.readlines():
        splitted = line[:-1].split(":")
        dict[splitted[0] + ":"] = stringComplex(splitted[1])
    f.close()
    return dict

def readMatrixChpr(m, path):
    s = matrixString(m)
    dict = readAllChpr(path)
    return dict[s]