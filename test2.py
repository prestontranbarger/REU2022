from SLTwoZ import *

q = 2 * 3

G = SL2Z
H = Gamma0(q)
inHChecker = inGroupChecker(H)
reps, gens, l, s = H.todd_coxeter()

#for rep in reps:
#    print(rep, ",")
#    print(TSDecomp(rep), ",")
#    print(contFrac(rep[0][0], rep[1][0]), "\n")

print(findCosetReps(reps, inHChecker, S * T ** 5 * S, q))
print(S * T ** 2 * S * T)