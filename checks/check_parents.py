"""Check orders, allowed signs and actual continued-fraction extensions exactly."""

from collections import Counter, defaultdict
from fractions import Fraction
from heapq import heappop, heappush
from itertools import takewhile
from math import lcm


def smooth_numbers(r, s):
    # Each row r^a, r^a s, ... is sorted. Open the next row at its predecessor.
    heap = [(1, 0)]
    while heap:
        value, b = heappop(heap)
        yield value
        heappush(heap, (value * s, b + 1))
        if b == 0:
            heappush(heap, (value * r, 0))


def unsigned_order(ell, v, base):
    if ell == 2:
        d = 1
        principal = base - (1 if base % 4 == 1 else -1)
    else:
        d = next(i for i in range(1, ell) if pow(base, i, ell) in (1, ell - 1))
        principal = base**d - (1 if pow(base, d, ell) == 1 else -1)
    c = 0
    while principal % ell == 0:
        principal //= ell
        c += 1
    return d * ell**max(0, v - c)


def group_generated(q, r, s):
    seen, todo = {1}, [1]
    while todo:
        x = todo.pop()
        for t in (x * r % q, x * s % q):
            if t not in seen:
                seen.add(t)
                todo.append(t)
    return seen


def expansions(p, q):
    digits = []
    x, y = p, q
    while y:
        d, rem = divmod(x, y)
        digits.append(d)
        x, y = y, rem
    assert digits[-1] >= 2
    return digits, digits[:-1] + [digits[-1] - 1, 1]


def convergents(digits):
    p0, q0, p1, q1 = 0, 1, 1, 0
    for digit in digits:
        p0, p1 = p1, digit * p1 + p0
        q0, q1 = q1, digit * q1 + q0
        yield p1, q1


def check_extension(p, q, successor, sign):
    parents = []
    matches = []
    for digits in expansions(p, q):
        conv = list(convergents(digits))
        assert conv[-1] == (p, q)
        h = conv[-2][1]
        parents.append(h)
        digit, remainder = divmod(successor - h, q)
        if remainder:
            continue
        assert digit >= 1
        # Interior point of the new cylinder even when the appended digit is 1.
        test = list(convergents(digits + [digit, 2]))
        assert test[-3] == (p, q) and test[-2][1] == successor
        next_num = test[-2][0]
        assert p * successor - q * next_num == sign
        ends = Fraction(next_num, successor), Fraction(next_num + p, successor + q)
        inside = Fraction(*test[-1])
        assert min(ends) < inside < max(ends)
        assert all(abs(t - Fraction(p, q)) <= Fraction(1, q * successor) for t in ends)
        matches.append(digits + [digit])
    assert sum(parents) == q and len(matches) == 1
    return matches[0]


def check_level(ell, v, r, s):
    q = ell**v
    group = group_generated(q, r, s)
    quotient = {min(x, q - x) for x in group}
    assert len(quotient) == lcm(unsigned_order(ell, v, r), unsigned_order(ell, v, s))
    if ell % 2:
        mod_ell = group_generated(ell, r, s)
        assert (q - 1 in group) == (ell - 1 in mod_ell)
    candidates = takewhile(lambda p: p < 2 * q, smooth_numbers(r, s))
    numerators = [p for p in candidates if p > q]
    targets = defaultdict(list)
    for p in numerators:
        for sign in (1, -1):
            residue = sign * pow(p, -1, q) % q
            if residue in group:
                targets[residue].append((p, sign))
    first_relation = next(m for m in smooth_numbers(r, s)
                          if m > 1 and m % q in (1, q - 1))
    counts = Counter()
    example = {}
    for m in smooth_numbers(r, s):
        if not targets:
            break
        if m > q:
            for p, sign in targets.pop(m % q, []):
                digits = check_extension(p, q, m, sign)
                assert p * m >= first_relation
                counts[sign] += 1
                if (ell, v, r, s, p) == (3, 4, 2, 5, 100):
                    example[sign] = (m, digits)
    assert counts[1] == len(numerators)
    assert counts[-1] == len(numerators) * (q - 1 in group)
    return counts, example


def check_nested_example():
    q, p = 3**10, 81920
    examined = 0
    for product in smooth_numbers(2, 5):
        examined += 1
        if product > 1 and product % q in (1, q - 1):
            break
    assert product == 2**59 * 5**8
    successor, remainder = divmod(product, p)
    assert remainder == 0 and successor == 2**45 * 5**7
    assert product % q == q - 1
    digits = check_extension(p, q, successor, -1)
    conv = list(convergents(digits))
    assert digits == [1, 2, 1, 1, 2, 1, 1, 4, 31, 1, 3, 2, 1, 46550814906941]
    assert conv[2:4] == [(4, 3), (7, 5)]
    assert conv[-2:] == [(p, q), ((product + 1) // q, successor)]
    assert Fraction(11, 8) < Fraction(p, q) < Fraction(7, 5)
    print('Nested example verified. Smooth products examined including 1:', examined)


if __name__ == '__main__':
    cases = [(2, 3, 5), (2, 17, 97), (3, 2, 5), (3, 19, 37),
             (5, 2, 11), (7, 2, 11), (7, 3, 5), (11, 3, 5)]
    total = Counter()
    levels = 0
    for ell, r, s in cases:
        v = 2 if ell == 2 else 1
        while ell**v <= 750:
            counts, example = check_level(ell, v, r, s)
            levels += 1
            total.update(counts)
            if example:
                print('q=81, p=100. Least smooth successor and expanded prefix by sign:', example)
            v += 1
    print('Levels:', levels, 'permitted triples by sign:', sorted(total.items()))
    check_nested_example()
