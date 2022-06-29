from SLTwoZ import *
import time

m = buildMatrix(7, 25)
print(m)

p = 5
G = SL2Z
H = Gamma0(p)
inHChecker = inGroupChecker(H)
reps = [matrix.identity(2),
        T ** 3 * S * T ** 4 * S * T ** 4,
        T ** 2 * S * T ** 4 * S,
        T ** 5 * S * T ** 5 * S * T * S * T ** 2 * S,
        T ** 6 * S * T ** 3 * S,
        T ** 4 * S * T ** 4 * S * T ** 4 * S * T ** 4 * S * T ** 4]
#reps = cosetReps(G, H, p)

bT = time.time()
tsD = TSDecomp(m)
reT = TSDecompToRewritingTape(tsD)
reW = reidemeisterSchreierRewriteReps(reps, inHChecker, reT, p)
tT = time.time() - bT

m = matrix.identity(2)
for letter in reW:
    print(letter, ",")
    m = m * letter
print("\nfinal:\n", m)
print(tT)