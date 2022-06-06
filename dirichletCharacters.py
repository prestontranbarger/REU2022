from base import *

def isDCharacter(dChar):
    #checks if a dictionary is a valid Dirichlet character
    ld = len(dChar)
    if dChar[0] != 0:
        return False
    for m in range(1, ld):
        d = dChar[m] + 0 * I
        if (float(d.real()) == 0 and float(d.imag()) == 0) != (gcd(m, ld) > 1):
            return False
        for n in range(ld):
            d1 = dChar[m * n % ld] + 0 * I
            d2 = dChar[m] * dChar[n] + 0 * I
            if (abs(float(d1.real()) - float(d2.real())) > 10 ** (-12)) or (abs(float(d1.imag()) - float(d2.imag())) > 10 ** (-12)):
                return False
    return True

def prodDCharacters(dChar1, dChar2):
    #multiplies two dirichlet characters
    if not dChar1:
        return dChar2
    ld1, ld2 = len(dChar1), len(dChar2)
    dChar = {}
    for i in range(lcm(ld1, ld2)):
        dChar[i] = dChar1[i % ld1] * dChar2[i % ld2]
    return dChar

def allDCharacters(n):
    #computes all dirichlet characters for a given modulus
    tn, fn = totient(n, True)
    if (n == 2 or n == 4) or (len(fn) == 1 and fn[0][0] != 2):
        pr = primitiveRoot(n)
        characters = []
        ztn = e ** (2 * pi * I / tn)
        for i in range(tn):
            character = {0: 0}
            for j in range(1, n):
                if gcd(j, n) > 1:
                    character[j] = 0
                character[pow(pr, j, n)] = ztn ** (i * j % tn)
            characters.append(character)
        return characters
    #elif len(fn) == 1 and fn[0][0] == 2:

    else:
        tCharacters = [allDCharacters(fn[i][0] ** fn[i][1]) for i in range(len(fn))]
        pCharacters = []
        for indxs in cartesianProd([[j for j in range(len(tCharacters[i]))] for i in range(len(tCharacters))]):
            character = {}
            for i in range(len(tCharacters)):
                character = prodDCharacters(character, tCharacters[i][indxs[i]])
            pCharacters.append(character)
        return pCharacters