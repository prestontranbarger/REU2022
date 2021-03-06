from sage.all import *
from tqdm import tqdm

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

def isPrime(n):
  return n in Primes()

def isPrimePower(n):
  f = factor(n)
  if len(f) == 1:
    return [True, f]
  return [False, f]

def pOppOc(n):
  [tf, f] = isPrimePower(n)
  if tf:
    if isPrime(n):
      return ['p', f[0][0]]
    else:
      return ['q', [f[0][0], f[0][1]]]
  else:
    return ['n']

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

def gammaFunction(z):
  return gamma(z)

def matrixString(m):
  return str(m[0][0]) + "," + str(m[0][1]) + "," + str(m[1][0]) + "," + str(m[1][1]) + ":"

def stringMatrix(s):
  l = s[:-1].split(",")
  return matrix(ZZ, [[int(l[0]), int(l[1])],
                     [int(l[2]), int(l[3])]])

def complexString(z):
  z += 0 * I
  return "(" + str(float(z.real())) + "," + str(float(z.imag())) + ")"

def stringComplex(s):
  l = s[1:-1].split(",")
  return float(l[0]) + float(l[1]) * I