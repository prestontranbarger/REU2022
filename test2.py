from SLTwoZ import *

p = 3
G = SL2Z
H = Gamma0(p)

l = []

reps = cosetReps(G, H, p)
inHChecker = inGroupChecker(H)

l.append(UReps(reps, inHChecker, findCosetReps(reps, inHChecker, matrix.identity(2), p), T, p))
l.append(UReps(reps, inHChecker, findCosetReps(reps, inHChecker, T, p), S, p))

m = matrix.identity(2)
for letter in l:
    m = m * letter
print(m)