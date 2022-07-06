from dirichletCharacters import *
from precomputation import *

for dChar1 in allDCharacters(4):
    for dChar2 in allDCharacters(7):
        if isEven(dChar1) == isEven(dChar2):
            if isPrimitive(dChar1) and isPrimitive(dChar2):
                print(dChar1, r'&', dChar2)
                precomputeCharacterPairs(dChar1, dChar2)