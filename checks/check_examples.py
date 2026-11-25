"""Reproduce the positive minima and the nested example in the paper."""

from fractions import Fraction


def least_signed_product(modulus, upper_bound):
    """Enumerate exactly the positive products 2^a 5^b up to upper_bound."""
    minimum = None
    count = 0
    power_of_two = 1
    while power_of_two <= upper_bound:
        product = power_of_two
        while product <= upper_bound:
            count += 1
            if product > 1 and product % modulus in (1, modulus - 1):
                if minimum is None or product < minimum:
                    minimum = product
            product *= 5
        power_of_two *= 2
    return minimum, count


def convergents(partial_quotients):
    previous_p, p = 0, 1
    previous_q, q = 1, 0
    result = []
    for term in partial_quotients:
        previous_p, p = p, term * p + previous_p
        previous_q, q = q, term * q + previous_q
        result.append((p, q))
    return result


def check():
    exponents = [(1, 0), (3, 0), (4, 1), (4, 1), (1, 7), (13, 10)]
    exponents.extend([(59, 8)] * 6)
    for level, (a, b) in enumerate(exponents, 1):
        candidate = 2**a * 5**b
        minimum, count = least_signed_product(3**level, candidate)
        assert minimum == candidate, (level, minimum, candidate)
        print(f"v={level:2d}  least N=2^{a} 5^{b}  ({count} products checked)")

    prefix = [1, 2, 1, 1, 2, 1, 1, 4, 31, 1, 3, 2, 1, 46550814906941]
    pairs = convergents(prefix)
    assert pairs[2:4] == [(4, 3), (7, 5)]
    numerator, denominator = 2**14 * 5, 3**10
    successor_q = 2**45 * 5**7
    product = numerator * successor_q
    assert product == 2**59 * 5**8
    assert product % denominator == denominator - 1
    successor_p, remainder = divmod(product + 1, denominator)
    assert remainder == 0
    assert pairs[-2:] == [(numerator, denominator), (successor_p, successor_q)]
    parent_p, parent_q = pairs[-3]
    assert parent_q == successor_q % denominator == 40891
    assert successor_p == parent_p + numerator * prefix[-1]
    assert abs(successor_p * denominator - product) == 1
    left, right = Fraction(11, 8), Fraction(7, 5)
    assert left < Fraction(numerator, denominator) < right
    for endpoint in (Fraction(successor_p, successor_q),
                     Fraction(successor_p + numerator, successor_q + denominator)):
        assert left < endpoint < right
    print(f"Nested final pair: {numerator}/{denominator}, {successor_p}/{successor_q}")


if __name__ == "__main__":
    check()
