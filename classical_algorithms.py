import time
import math
import random
from math import gcd

def pollards_rho(n):
    def f(x):
        return (x**2 + 1) % n 
    x = random.randint(2, n - 1)
    y = x
    d = 1
    while d == 1:
        x = f(x)
        y = f(f(y))
        d = gcd(abs(x - y), n)
    return d if d != n else None

def factors_brute_force(n):
    factors = []
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)
    return factors

def measure_time(function, n):
    start_time = time.time()
    result = function(n)
    elapsed_time = time.time() - start_time
    return result, elapsed_time

if __name__ == "__main__":
    number = 314191  

    print("\n--- Factors using different algorithms ---")
    factors, time_taken = measure_time(factors_brute_force, number)
    print(f"Brute force: {factors}")
    print(f"Time taken: {time_taken:.6f} seconds\n")
    factors, elapsed_time = measure_time(pollards_rho, number)
    print(f"Pollard's Rho of {number}: {factors}")
    print(f"Time taken for factorization: {elapsed_time:.6f} seconds")
