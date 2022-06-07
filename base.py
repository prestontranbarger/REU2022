from sage.all import *

def arrToDict(arr):
  dict = {}
  for i in range(len(arr)):
    dict[i] = arr[i]
  return dict

def cartesianProd(pools):
  #computes the cartesian product of arrays
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

def jacobiSymbol(a, n):
  #generalized Legendre symbol
  if n == 1:
    return [1, True]
  if n % 2 == 0 or n < 0:
    return [0, False]
  if math.gcd(a, n) != 1:
    return [0, True]
  else:
    prod = 1
    if a < 0:
      a *= -1
      if n % 4 == 3:
        prod *= -1
    while a > 1:
      a %= n
      m = (-1 if (n % 8 == 3 or n % 8 == 5) else 1)
      while not a % 2:
        a = a // 2
        prod *= m
      if a > 2:
        if (a - 1) * (n - 1) / 4 % 2:
          prod *= -1
        a, n = n, a
    return [prod, True]