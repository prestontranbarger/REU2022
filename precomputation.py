from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

import os
from os.path import isfile
from os.path import isdir

gensWritePath = "/precomputation"
charWritePath = "/precomputation"

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

def precomputeGensGammaZero(n):
    subPath = os.getcwd() + gensWritePath + str("/gammaZeroGens/")
    filePath = subPath + "G0-" + str(n) + ".gens"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        gens = [stringMatrix(line[:-1]) for line in f.readlines()]
        f.close()
    else:
        reps = cosetRepsSLTwoZOverGammaZeroGeneralized(n)
        inHChecker = inGroupChecker(Gamma0(n))
        gens = schreierLemmaReps(SL2Z, reps, inHChecker, n)
        f = open(filePath, "w")
        f.writelines([matrixString(gen) + "\n" for gen in gens])
        f.close()
    return gens

def precomputeGensGammaOne(n):
    subPath = os.getcwd() + gensWritePath + str("/gammaOneGens/")
    filePath = subPath + "G1-" + str(n) + ".gens"
    if not isdir(subPath):
        os.mkdir(subPath)
    if os.path.exists(filePath):
        f = open(filePath, "r")
        gens = [stringMatrix(line[:-1]) for line in f.readlines()]
        f.close()
    else:
        reps = cosetRepsSLTwoZOverGammaOneGeneralized(n)
        inHChecker = inGroupChecker(Gamma1(n))
        gens = schreierLemmaReps(SL2Z, reps, inHChecker, n)
        f = open(filePath, "w")
        f.writelines([matrixString(gen) + "\n" for gen in gens])
        f.close()
    return gens

def precomputeCharacterPairs(dChar1, dChar2):
    q1q2 = modulus(dChar1) * modulus(dChar2)
    subPath, (G0Path, G1Path) = createCharacterPairFiles(dChar1, dChar2)
    f = open(G0Path, "w")
    f.writelines([matrixString(gen) + complexString(newFormDedekindSum(dChar1, dChar2, gen)) + "\n" for gen in tqdm(precomputeGensGammaZero(q1q2))])
    f.close()
    f = open(G1Path, "w")
    f.writelines([matrixString(gen) + complexString(newFormDedekindSum(dChar1, dChar2, gen))+ "\n" for gen in tqdm(precomputeGensGammaOne(q1q2))])
    f.close()

def createCharacterPairFiles(dChar1, dChar2):
    q1q2 = modulus(dChar1) * modulus(dChar2)
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
    G0Path,\
    G1Path = subPath + "G0-" + str(q1q2) + ".chpr",\
             subPath + "G1-" + str(q1q2) + ".chpr"
    f = open(G0Path, "w")
    f.close()
    f = open(G1Path, "w")
    f.close()
    return subPath, (G0Path, G1Path)