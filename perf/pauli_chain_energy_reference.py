"""Independent numerical Heisenberg-chain reference; requires NumPy/SciPy.

Even isotropic spin-1/2 chains have their ground energy in total Sz=0:
every integer-total-spin multiplet contains a zero-magnetization member.
This is a numerical eigenvalue reference, not an interval certificate.
"""
import json
import time

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def ground_energy(n):
    started = time.perf_counter()
    states = np.arange(1 << n, dtype=np.uint32)
    states = states[np.bitwise_count(states) == n // 2]
    dim = len(states)
    index = np.full(1 << n, -1, dtype=np.int32)
    index[states] = np.arange(dim, dtype=np.int32)
    rows, cols, values = [], [], []
    diagonal = np.zeros(dim)
    for i in range(n):
        j = (i + 1) % n
        opposite = ((states >> i) & 1) != ((states >> j) & 1)
        diagonal += np.where(opposite, -0.25, 0.25)
        selected = np.flatnonzero(opposite)
        rows.append(selected)
        cols.append(index[states[selected] ^ ((1 << i) | (1 << j))])
        values.append(np.full(len(selected), 0.5))
    rows.append(np.arange(dim))
    cols.append(np.arange(dim))
    values.append(diagonal)
    h = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(dim, dim)).tocsr()
    assert (h - h.T).nnz == 0
    v0 = np.random.default_rng(42).normal(size=dim)
    if dim <= 4:
        energies, vectors = np.linalg.eigh(h.toarray())
        energy, vector = energies[0], vectors[:, 0]
    else:
        energies, vectors = eigsh(h, k=1, which="SA", v0=v0, tol=1e-12)
        energy, vector = energies[0], vectors[:, 0]
    residual = np.linalg.norm(h @ vector - energy * vector)
    assert residual <= 1e-9
    if n in (2, 4):
        assert abs(energy - {2: -1.5, 4: -2.0}[n]) <= 1e-12
    return dict(n=n, dimension=dim, total_energy=float(energy), per_spin=float(energy/n),
                residual=float(residual), seconds=time.perf_counter()-started)


if __name__ == "__main__":
    print(json.dumps(dict(numpy=np.__version__, scipy=scipy.__version__)), flush=True)
    for size in (2, 4, 10, 20):
        print(json.dumps(ground_energy(size)), flush=True)
