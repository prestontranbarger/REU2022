from SLTwoZ import *

m = matrix(ZZ, [[10, 11],
                [ 9, 10]])
print(TSDecomp(m), "\n")
n = 9
reps = cosetRepsSLTwoZOverGammaOneGeneralized(n)
inHChecker = inGroupChecker(Gamma1(n))
t0 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, matrix.identity(2), n), -1 * matrix.identity(2), n)
t1 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, -1 * matrix.identity(2), n), T, n)
t2 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, -1 * T, n), S, n)
t3 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, -1 * T * S * T ** (-9), n), T ** (-9), n)
t4 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, -1 * T * S * T ** (-9), n), S, n)
t5 = UReps(reps, inHChecker, findCosetReps(reps, inHChecker, -1 * T * S * T ** (-9) * S, n), T, n)
print(t0)
print(t1)
print(t2)
print(t3)
print(t4)
print(t5)
print(t0 * t1 * t2 * t3 * t4 * t5)