"""Independent exact-arithmetic checks for The size of a smooth neighboring convergent."""
from collections import deque
from fractions import Fraction
from math import gcd, isqrt


def cf_digits(p, q):
    digits = []
    while q:
        digit, remainder = divmod(p, q)
        digits.append(digit)
        p, q = q, remainder
    return digits


def convergents(digits):
    a, b = 0, 1
    c, d = 1, 0
    result = []
    for digit in digits:
        a, c = c, a + digit * c
        b, d = d, b + digit * d
        result.append((c, d))
    return result


def check_cylinders(max_q=60):
    # Start from all rational successors in the range, independent of the congruence.
    observed = set()
    for Q in range(4, 3 * max_q + 1):
        for successor in range(Q + 1, 2 * Q):
            if gcd(successor, Q) != 1:
                continue
            canonical = cf_digits(successor, Q)
            split = canonical[:-1] + [canonical[-1] - 1, 1]
            for digits in (canonical, split):
                (P, q), final = convergents(digits)[-2:]
                assert final == (successor, Q)
                if not (3 <= q <= max_q and q < P < 2 * q and q < Q <= 3 * q):
                    continue
                observed.add((P, q, Q))
                x = Fraction(P, q)
                near = Fraction(successor + P, Q + q)
                far = Fraction(successor, Q)
                assert abs(x - near) == Fraction(1, q * (Q + q))
                assert abs(x - far) == Fraction(1, q * Q)
                assert (near - x) * (far - x) > 0
    tested = 0
    for q in range(3, max_q + 1):
        for P in range(q + 1, 2 * q):
            if gcd(P, q) != 1:
                continue
            for Q in range(q + 1, 3 * q + 1):
                tested += 1
                assert ((P, q, Q) in observed) == (P * Q % q in (1, q - 1))
    print(f"Cylinder check: {tested} reduced-numerator triples, {len(observed)} successes.")


def valuation(n, ell):
    assert n != 0
    n = abs(n)
    count = 0
    while n % ell == 0:
        n //= ell
        count += 1
    return count


def check_valuations():
    trials = 0
    for ell in (2, 3, 5, 7, 11):
        h = 2 * (ell - 1)
        for a, b in [(2, 3), (3, 5), (5, 3), (11, 7), (9, 17), (19, 1)]:
            if (a * b) % ell == 0:
                continue
            e = valuation(a**h - b**h, ell)
            for m in range(1, 45):
                assert valuation(a**(h*m) - b**(h*m), ell) == e + valuation(m, ell)
                trials += 1
    print(f"Unit valuations: {trials} exact rational-unit and exponent tests.")


def signed_key(value, q):
    t = value % q
    return min(t, q - t)


def check_mixed_moduli():
    tested = 0
    for p, s, r in [(2, 3, 5), (2, 5, 3), (3, 5, 2), (3, 7, 11)]:
        Cp = p**valuation(r**(2*(p-1))-1, p)
        Cs = s**valuation(r**(2*(s-1))-1, s)
        for i in range(1, 7):
            for j in range(1, 6):
                q = p**i * s**j
                if q > 18000:
                    continue
                value = r % q
                m = 1
                while signed_key(value, q) != 1:
                    value = r*value % q
                    m += 1
                assert (m * Cp * Cs) % q == 0
                tested += 1
    print(f"Mixed-modulus corollary: first signed orders checked for {tested} moduli and bases.")


def signed_index(q, r, s):
    seen = {1}
    queue = deque([1])
    while queue:
        t = queue.popleft()
        for prime in (r, s):
            k = signed_key(t * prime, q)
            if k not in seen:
                seen.add(k)
                queue.append(k)
    return len(seen)


def ceil_sqrt(n):
    return isqrt(n) + (isqrt(n)**2 != n)


def shortest(q, r, s, index):
    # The pigeonhole square gives a certified Euclidean search radius.
    radius = ceil_sqrt(2 * ceil_sqrt(index)**2)
    powers_s = {b: pow(s, b, q) for b in range(-radius, radius + 1)}
    candidates = []
    for a in range(radius + 1):
        ra = pow(r, a, q)
        for b, sb in powers_s.items():
            if a == 0 and b <= 0:
                continue
            if signed_key(ra * sb, q) == 1:
                candidates.append((a*a + b*b, a, b))
    return min(candidates)[1:]


def det(u, v):
    return u[0]*v[1] - u[1]*v[0]


def smooth_pairs(r, s, limit):
    power_s = 1
    b = 0
    while power_s <= limit:
        value = power_s
        a = 0
        while value <= limit:
            yield value, a, b
            value *= r
            a += 1
        power_s *= s
        b += 1


def check_lattices():
    levels = changes = constructions = 0
    for ell, r, s in [(2, 3, 5), (2, 5, 7), (3, 2, 5), (3, 7, 11),
                      (5, 2, 3), (7, 2, 11), (11, 3, 7)]:
        q = ell
        prev = None
        while q <= 2200:
            if q == 2:
                q *= ell
                continue
            index = signed_index(q, r, s)
            z = shortest(q, r, s, index)
            levels += 1
            if prev is not None:
                y = (ell * prev[0], ell * prev[1])
                assert signed_key(pow(r, y[0], q)*pow(s, y[1], q), q) == 1
                determinant = det(z, y)
                assert determinant % index == 0
                if determinant:
                    changes += 1
                    # Build a positive solution near a prescribed numerator by exact rounding.
                    radius = (ceil_sqrt(sum(t*t for t in z)) + ceil_sqrt(sum(t*t for t in y)) + 1)//2
                    for P, u, w in smooth_pairs(r, s, 3*q//2):
                        if not (4*q < 3*P and 2*P < 3*q):
                            continue
                        target = (u + 2*radius, w + 2*radius)
                        m = round(Fraction(det(target, y), determinant))
                        n = round(Fraction(det(z, target), determinant))
                        A, B = (m*z[0] + n*y[0], m*z[1] + n*y[1])
                        assert (A-target[0])**2 + (B-target[1])**2 <= radius**2
                        assert A-u >= radius and B-w >= radius
                        Q = r**(A-u)*s**(B-w)
                        assert Q > q and signed_key(P*Q, q) == 1
                        constructions += 1
                        break
            prev = z
            q *= ell
    print(f"Seven prime triples: {levels} levels, {changes} direction changes, "
          f"{constructions} explicit positive quadrant constructions.")


def check_three_adic():
    # 4^theta * 25 = 1 mod 3^v. Lift theta modulo 3^(v-1) one digit at a time.
    theta = 0
    tests = 0
    for v in range(2, 13):
        q = 3**v
        possible = [theta + j*(q//9) for j in range(3)]
        (theta,) = [t for t in possible if pow(4, t, q)*25 % q == 1]
        for a in range(-12, 13):
            for b in range(-12, 13):
                lhs = signed_key(pow(2, a, q) * pow(5, b, q), q) == 1
                assert lhs == ((a - b*theta) % (q//3) == 0)
                tests += 1
    print(f"3-adic ratio: {tests} signed exponent comparisons at levels 2 through 12.")


def least_product(q):
    # A residue table supplies the least possible exponent of 2 for each exponent of 5.
    logs = {}
    residue = 1
    while residue not in logs:
        logs[residue] = len(logs)
        residue = 2*residue % q
    period = len(logs)
    best = 2**logs[q-1]
    five = 1
    while five <= best:
        for sign in (1, -1):
            target = sign*pow(five, -1, q) % q
            if target not in logs:
                continue
            a = logs[target]
            if five == 1 and a == 0:
                a = period
            best = min(best, (1 << a)*five)
        five *= 5
    return best, logs[q-1], period


def interval_extremum(q, Q0, lo, hi):
    # Any contender with Q <= Q0 has P*Q < hi*q*Q0.
    cap = (hi*q*Q0).__ceil__() - 1
    pairs = []
    scanned = 0
    signed = 0
    for n, a, b in smooth_pairs(2, 5, cap):
        scanned += 1
        if n <= 1 or signed_key(n, q) != 1:
            continue
        signed += 1
        for i in range(a + 1):
            for j in range(b + 1):
                P = 2**i * 5**j
                Q = n//P
                if not (lo*q < P < hi*q and q < Q <= Q0):
                    continue
                successor = (n+1)//q if n%q == q-1 else (n-1)//q
                left, right = sorted([Fraction(successor, Q), Fraction(successor+P, Q+q)])
                if max(left, lo) < min(right, hi):
                    pairs.append((Q, P))
    return min(pairs), scanned, signed


def check_example():
    q = 3**10
    N, minus_exponent, period = least_product(q)
    assert (minus_exponent, period) == (19683, 39366)
    assert N == 2**59*5**8
    size = sum(1 for _ in smooth_pairs(2, 5, N))
    residues = sorted(n for n, a, b in smooth_pairs(2, 5, N) if signed_key(n,q)==1)
    assert size == 1352 and residues == [1, N]
    P = 81920
    Q = N//P
    successor = (N+1)//q
    assert Q == 2**45*5**7 and successor*q - P*Q == 1
    digits = cf_digits(successor, Q)
    assert digits == [1,2,1,1,2,1,1,4,31,1,3,2,1,46550814906941]
    cs = convergents(digits)
    assert cs[2:4] == [(4,3),(7,5)] and cs[-2:] == [(P,q),(successor,Q)]
    assert Fraction(11,8) < Fraction(P,q) < Fraction(7,5)
    print(f"Exact example: q={q}, N={N}, P={P}, Q={Q}, P'={successor}.")
    print(f"Order of 2: {period}. Products through N including 1: {size}. "
          f"Only nontrivial signed solution: {N}.")
    print("Successor Euclidean expansion:", digits)
    for lo, hi, expected in [(Fraction(1), Fraction(2), (2**47*5**6, 102400)),
                             (Fraction(11,8), Fraction(7,5), (Q, P))]:
        pair, count, signed = interval_extremum(q, Q, lo, hi)
        assert pair == expected
        print(f"Interval ({lo},{hi}): least Q={pair[0]}, P={pair[1]}. "
              f"Certified using {count} products, {signed} nontrivial signed solutions.")
    # Compare both supports at exactly the same rational numerator and denominator.
    fixed = 2**16
    q_single = 2**(minus_exponent-16)
    q_two = N//fixed
    assert q < fixed < 2*q and q < q_two and N%fixed == 0
    assert 10**5920 <= q_single < 10**5921 and 10**18 <= q_two < 10**19
    print(f"Fixed P=65536: one prime gives Q=2^19667 (5921 digits); "
          f"two give Q=2^43*5^8={q_two} (19 digits).")


if __name__ == '__main__':
    check_cylinders()
    check_valuations()
    check_mixed_moduli()
    check_lattices()
    check_three_adic()
    check_example()

