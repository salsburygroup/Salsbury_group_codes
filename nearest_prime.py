import argparse
from sympy import factorint

def nearest_prime_product(n, primes):
    # Start with the given number and search outwards
    lower = upper = n

    while True:
        # Check the lower number
        factors = factorint(lower)
        if all(factor in primes for factor in factors):
            return lower

        # Check the upper number
        factors = factorint(upper)
        if all(factor in primes for factor in factors):
            return upper

        # Expand the search
        lower -= 1
        upper += 1

def main():
    parser = argparse.ArgumentParser(description='Find the nearest integer that can be factored into powers and products of given prime numbers.')
    parser.add_argument('n', type=int, help='The input integer')
    parser.add_argument('primes', type=int, nargs='+', help='A set of prime numbers')
    args = parser.parse_args()

    result = nearest_prime_product(args.n, set(args.primes))
    print(result)

if __name__ == "__main__":
    main()

