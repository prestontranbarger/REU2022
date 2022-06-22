from SLTwoZ import *
import time

m = T * S * T ** 2 * S * T ** 4 * S * T ** 3 * S
print(m)

p = 19
G = SL2Z
H = Gamma0(p)
inHChecker = inGroupChecker(H)

bT = time.time()
tsD = TSDecomp(m)
reT = TSDecompToRewritingTape(tsD)
reps = cosetReps(G, H, p)
reW = reidemeisterSchreierRewriteReps(reps, inHChecker, reT, p)

m = matrix.identity(2)
for letter in reW:
    print(letter, ",")
    m = m * letter
print("\nfinal:\n", m)
print(float(time.time() - bT))