from sage.all import *

def cartesianProd(pools):
  result = [[]]
  for pool in pools:
    result = [x + [y] for x in result for y in pool]
  return result

def totient(n, rf = False):
  # computes euler totient function of n, rf is an optional parameter to return the factorization of n
  f = factor(n)
  p = 1
  for pe in f:
    p *= (pe[0] - 1) * (pe[0] ** (pe[1] - 1))
  if rf:
    return p, f
  else:
    return p


def primitiveRoot(n):
  # computes primitive roots if they exist otherwise, returns 0
  tn = totient(n)
  ftn = factor(tn)
  for g in range(1, n):
    if gcd(n, g) == 1:
      flag = True
      for i in range(len(ftn)):
        if pow(g, tn // ftn[i][0], n) == 1:
          flag = False
      if flag:
        return g
  return 0