# -*- coding: utf-8 -*-

"""Tests of the pure sampling helpers on a synthetic water box (no molsystem)."""

import numpy as np
import pytest

from extract_clusters_step.cluster_sampling import (
    classify_motif,
    grow_connected,
    molecule_neighbour_graph,
    unwrap_shifts,
)

L = 23.0  # cubic box, Å (~1 g/cc for 400 waters)
N = 400
T = np.eye(3) * L
CELL = (L, L, L, 90.0, 90.0, 90.0)


@pytest.fixture(scope="module")
def water_box():
    """400 rigid waters at random positions; whole molecules may exceed [0,1)."""
    rng = np.random.default_rng(1)
    oxy = rng.random((N, 3))  # fractional
    h1 = np.array([0.7572, 0.5865, 0.0]) / L
    h2 = np.array([-0.7572, 0.5865, 0.0]) / L
    frac = np.concatenate([np.stack([o, o + h1, o + h2]) for o in oxy])
    atom_mol = np.repeat(np.arange(N), 3)
    symbols = np.array(["O", "H", "H"] * N)
    return dict(O=oxy, frac=frac, atom_mol=atom_mol, mask_O=symbols == "O", rng=rng)


@pytest.fixture(scope="module")
def adj(water_box):
    b = water_box
    return molecule_neighbour_graph(
        b["frac"], b["atom_mol"], b["mask_O"], T, 3.5, True, CELL
    )


def _mic(d):
    d = d - np.round(d)
    return np.linalg.norm(d @ T)


def test_graph_paths_agree(water_box, adj):
    """The orthorhombic KD-tree path and the general MIC path give the same graph."""
    b = water_box
    # force the general path on the same geometry by lying about an angle
    adj_gen = molecule_neighbour_graph(
        b["frac"],
        b["atom_mol"],
        b["mask_O"],
        T,
        3.5,
        True,
        (L, L, L, 90.0, 90.0, 90.0 + 1e-3),
    )
    assert adj == adj_gen
    n_edges = sum(len(v) for v in adj.values()) // 2
    assert n_edges > 0


def test_edges_match_brute_force(water_box, adj):
    """Edges are exactly the O-O pairs within the cutoff under minimum image."""
    oxy = water_box["O"]
    for i in list(adj)[:50]:
        for j in adj[i]:
            assert _mic(oxy[j] - oxy[i]) < 3.5 + 1e-9
        for j in range(N):
            if j != i and j not in adj[i]:
                assert _mic(oxy[j] - oxy[i]) >= 3.5 - 1e-9


@pytest.mark.parametrize("n", [3, 4, 5])
def test_growth_connectivity_and_unwrap(water_box, adj, n):
    """Grown clusters are connected and contiguous after unwrapping."""
    rng = np.random.default_rng(n)
    oxy = water_box["O"]
    n_ok = 0
    for _ in range(300):
        g = grow_connected(adj, int(rng.integers(N)), n, rng)
        if g is None:
            continue
        order, parent = g
        assert len(set(order)) == n
        # connected in the induced subgraph
        inside = set(order)
        seen = {order[0]}
        stack = [order[0]]
        while stack:
            u = stack.pop()
            for v in adj[u] & inside:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        assert seen == inside
        # after unwrapping every internal contact pair is contiguous
        shifts = unwrap_shifts(order, parent, oxy, True)
        xyz = {m: (oxy[m] + shifts[m]) @ T for m in order}
        for i in order:
            for j in adj[i] & inside:
                d_plain = np.linalg.norm(xyz[i] - xyz[j])
                assert abs(d_plain - _mic(oxy[i] - oxy[j])) < 1e-9
        n_ok += 1
    assert n_ok > 100


def test_grow_returns_none_for_small_component():
    adj = {0: {1}, 1: {0}}
    assert grow_connected(adj, 0, 3, np.random.default_rng(0)) is None
    assert grow_connected({}, 5, 2, np.random.default_rng(0)) is None


def test_classify_motif():
    assert classify_motif(2, 1, [1, 1]) == "dimer"
    assert classify_motif(3, 2, [2, 1, 1]) == "chain"
    assert classify_motif(3, 3, [2, 2, 2]) == "ring"
    assert classify_motif(4, 3, [3, 1, 1, 1]) == "star"
    assert classify_motif(4, 3, [2, 2, 1, 1]) == "chain"
    assert classify_motif(4, 4, [2, 2, 2, 2]) == "ring"
    assert classify_motif(4, 4, [3, 2, 2, 1]) == "paw"
    assert classify_motif(4, 5, [3, 3, 2, 2]) == "diamond"
    assert classify_motif(4, 6, [3, 3, 3, 3]) == "K4"
    assert classify_motif(6, 7, [3, 3, 2, 2, 2, 2]) == "e7"


def test_motif_names():
    from extract_clusters_step.cluster_sampling import motif_names

    assert motif_names(2) == ["dimer"]
    assert motif_names(3) == ["chain", "ring"]
    assert motif_names(4) == ["chain", "star", "ring", "paw", "diamond", "K4"]
    assert motif_names(5) == ["e4", "e5", "e6", "e7", "e8", "e9", "e10"]
    # every label classify_motif can produce for n=4 is in the list
    for n_edges, degrees in (
        (3, [3, 1, 1, 1]),
        (3, [2, 2, 1, 1]),
        (4, [2, 2, 2, 2]),
        (4, [3, 2, 2, 1]),
        (5, [3, 3, 2, 2]),
        (6, [3, 3, 3, 3]),
    ):
        assert classify_motif(4, n_edges, degrees) in motif_names(4)
