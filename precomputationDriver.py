from precomputation import *

for n in tqdm(range(3, 31)):
    if n not in Primes():
        precomputeGensGammaZero(n)
        precomputeGensGammaOne(n)