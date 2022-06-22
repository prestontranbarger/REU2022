from SLTwoZ import *
import time

m = matrix(ZZ, [[ 8, -3],
                [19, -7]])
m = matrix(ZZ, [[      52,      19],
                [23456773, 8570744]])
print(m)

p = 19
G = SL2Z
H = Gamma0(p)
inHChecker = inGroupChecker(H)
reps = cosetReps(G, H, p)

bT = time.time()
tsD = TSDecomp(m)
tT = time.time() - bT
reT = TSDecompToRewritingTape(tsD)
reW = reidemeisterSchreierRewriteReps(reps, inHChecker, reT, p)

m = matrix.identity(2)
for letter in reW:
    print(letter, ",")
    m = m * letter
print("\nfinal:\n", m)
print(tT)