from additiveCharacters import *
from dedekindSum import *
from dirichletCharacters import *
from SLTwoZ import *

dChar1 = allDCharacters(5)[3]
print(getValues(dChar1))
dChar2 = allDCharacters(7)[5]
print(getValues(dChar2))
print(getValues(prodDCharacters(dChar1, dChar2)))
