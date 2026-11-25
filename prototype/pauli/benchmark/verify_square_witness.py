"""Strict exact ED in every Mz sector; requires python-flint (not NumPy)."""
import re
import sys
from collections import defaultdict
from itertools import product

from flint import fmpz_mat


def verify(path):
    with open(path) as f:
        scale, floor, denominator = map(int, re.fullmatch(
            r"# scale=(\d+) floor=(-?\d+) denominator=(\d+)\n?", f.readline()).groups())
        if scale != denominator or scale.bit_count() != 1:
            raise ValueError("Expected a common power-of-two denominator")
        terms = {}
        for line in f:
            *key, v = map(int, line.split())
            if len(key) == 4 and key[2:] == [0, 0]:
                key = key[:2]
            pairs = list(zip(key[::2], key[1::2]))
            if (len(key) not in (2, 4, 6) or len(set(key)) != len(key)
                    or any(not 1 <= i < j <= 9 for i, j in pairs)
                    or pairs != sorted(pairs) or tuple(key) in terms):
                raise ValueError(f"Invalid or repeated observable: {key}")
            terms[tuple(key)] = v
    if not terms:
        raise ValueError("Empty cut")

    # Build exact Pauli X/Y/Z actions, independently of the swap expansion
    # used to propose the floor. Product phases are real (every Y is paired).
    masks = defaultdict(int)
    for sites, coef in terms.items():
        for pair_axes in product(range(3), repeat=len(sites)//2):
            flip = sign = ny = 0
            for site, a in zip(sites, (a for a in pair_axes for _ in range(2))):
                bit = 1 << (site-1)
                flip |= bit if a != 2 else 0
                sign |= bit if a != 0 else 0
                ny += a == 1
            masks[(flip, sign)] += coef*(-1)**(ny//2)
    entries = [[0]*512 for _ in range(512)]
    for (flip, sign), coef in masks.items():
        if coef:
            for s in range(512):
                entries[s ^ flip][s] += coef*(-1)**((s & sign).bit_count() & 1)
    for i in range(512):
        for j in range(i):
            if entries[i][j] != entries[j][i] or (
                    i.bit_count() != j.bit_count() and entries[i][j]):
                raise ArithmeticError("The witness is non-symmetric or mixes Mz sectors")
    for weight in range(10):
        states = [i for i in range(512) if i.bit_count() == weight]
        n = len(states)
        # det(tI + scale*H - floor*I) has real roots. Strictly positive
        # coefficients preclude any nonnegative root, proving positive definiteness.
        neg = fmpz_mat(n, n, [floor*(i == j)-entries[i][j] for i in states for j in states])
        if not all(c > 0 for c in neg.charpoly()):
            raise ArithmeticError(f"Floor failed exact positive definiteness in sector {weight}")
        print(f"Mz weight {weight}: positive definite, dimension {n}")
    print(f"Exact floor {floor}/{scale}; {len(terms)} nonzero observable coefficients")


if __name__ == "__main__":
    verify(sys.argv[1])
