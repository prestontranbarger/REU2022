from base import *

def isACharacter(aChar):
    #checks if a dictionary is a valid additive character
    ld = len(aChar)
    for m in range(ld):
        for n in range(ld):
            d1 = aChar[(m + n) % ld] + 0 * I
            d2 = aChar[m] * aChar[n] + 0 * I
            if (abs(float(d1.real()) - float(d2.real())) > 10 ** (-12)) or (abs(float(d1.imag()) - float(d2.imag())) > 10 ** (-12)):
                return False
    return True

def prodACharacters(aChar1, aChar2):
    #multiplies two additive characters
    if not aChar1:
        return aChar2
    ld1, ld2 = len(aChar1), len(aChar2)
    aChar = {}
    for i in range(ld1 * ld2):
        aChar[i] = aChar1[i % ld1] * aChar2[i % ld2]
    return aChar

def allACharacters(n):
    #computes all additive characters for a given modulus
    characters = []
    zn = e ** (2 * pi * I / n)
    for i in range(n):
        character = {}
        for j in range(n):
            character[j] = zn ** (i * j % n)
        characters.append(character)
    return characters