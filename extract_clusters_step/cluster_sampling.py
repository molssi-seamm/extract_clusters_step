# -*- coding: utf-8 -*-

"""Sampling of n-molecule clusters from a (periodic) configuration.

The pure functions here have no dependence on the SEAMM framework, only on
numpy/scipy and the duck-typed ``molsystem`` configuration passed in, so they
can be unit-tested directly and reused outside a flowchart.

The approach:

1. Find the molecules (by bonds) and build a molecular **contact graph**:
   molecule A is adjacent to B if any *contact atom* of A is within ``cutoff``
   (minimum image) of any contact atom of B.
2. Sample n-mers as **connected subgraphs** of that graph by random frontier
   growth from a random seed molecule. Every member interacts with at least one
   other, and chains, rings and stars all occur -- unlike nearest-neighbour
   selection, which only ever yields the most compact cluster.
3. Optionally **stratify** the accepted clusters to be flat in a spread
   coordinate (radius of gyration or largest centroid separation of the
   molecules), with bin edges taken as quantiles of a pilot sample, and
   optionally balance jointly over contact-graph topology (motif).
4. Emit each cluster as an **unwrapped, centred, non-periodic** configuration
   with the molecules intact, a unique name carrying provenance, and optional
   per-configuration properties.
"""

import math
from collections import Counter
import logging

import numpy as np

logger = logging.getLogger(__name__)

PROPERTY_TAG = "#ExtractClusters#scan"

SPREAD_METRICS = ("rg", "dmax")


def classify_motif(n, n_edges, degrees):
    """Name the topology of the induced contact graph of an n-mer.

    Only n <= 4 get descriptive names; larger clusters are labelled by edge
    count so they can still be balanced on topology.

    Parameters
    ----------
    n : int
        Number of molecules.
    n_edges : int
        Number of contact edges within the cluster.
    degrees : [int]
        The degree of each molecule within the cluster.

    Returns
    -------
    str
    """
    maxdeg = max(degrees) if degrees else 0
    if n == 2:
        return "dimer"
    if n == 3:
        return {2: "chain", 3: "ring"}.get(n_edges, f"e{n_edges}")
    if n == 4:
        if n_edges == 3:
            return "star" if maxdeg == 3 else "chain"
        if n_edges == 4:
            return "ring" if maxdeg == 2 else "paw"  # triangle + tail
        if n_edges == 5:
            return "diamond"
        if n_edges == 6:
            return "K4"
    return f"e{n_edges}"


def molecule_neighbour_graph(
    frac_atoms, atom_mol, contact_mask, T, cutoff, periodic, cell_parameters
):
    """Molecule-level adjacency from atom contacts.

    A~B if any contact atom of A is within ``cutoff`` (Å, minimum image) of any
    contact atom of B.

    Parameters
    ----------
    frac_atoms : (N,3) ndarray
        Fractional coordinates (Cartesian if not periodic).
    atom_mol : (N,) ndarray of int
        Molecule index of each atom.
    contact_mask : (N,) ndarray of bool
        Which atoms define contact.
    T : (3,3) ndarray
        Fractional->Cartesian transform (``xyz = uvw @ T``); identity if not
        periodic.
    cutoff : float
        Contact distance, Å.
    periodic : bool
    cell_parameters : (a, b, c, alpha, beta, gamma) or None

    Returns
    -------
    dict[int, set[int]]
        Adjacency lists keyed by molecule index. Molecules with no contacts
        are absent.
    """
    from scipy.spatial import cKDTree

    idx = np.nonzero(contact_mask)[0]
    mol_of = atom_mol[idx]

    if not periodic:
        xyz = frac_atoms[idx]
        pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
    else:
        a, b, c, alpha, beta, gamma = cell_parameters
        orthorhombic = all(abs(x - 90.0) < 1e-6 for x in (alpha, beta, gamma))
        if orthorhombic:
            # Fast path: periodic KD-tree on wrapped Cartesian coordinates.
            L = np.array([a, b, c])
            xyz = (frac_atoms[idx] % 1.0) * L
            xyz = np.minimum(xyz, np.nextafter(L, 0))  # guard x == L
            pairs = cKDTree(xyz, boxsize=L).query_pairs(cutoff, output_type="ndarray")
        else:
            # General cell: exact minimum-image on fractional differences,
            # chunked so memory stays bounded. O(M^2) in contact atoms.
            f = frac_atoms[idx]
            M = len(f)
            out = []
            chunk = max(1, int(2e7 // max(M, 1)))
            for s in range(0, M, chunk):
                d = f[s : s + chunk, None, :] - f[None, :, :]
                d -= np.round(d)
                r = np.linalg.norm(d @ T, axis=2)
                ii, jj = np.nonzero(r < cutoff)
                ii += s
                keep = ii < jj
                out.append(np.stack([ii[keep], jj[keep]], 1))
            pairs = np.concatenate(out) if out else np.zeros((0, 2), int)

    adj = {}
    if len(pairs):
        mi, mj = mol_of[pairs[:, 0]], mol_of[pairs[:, 1]]
        keep = mi != mj
        for i, j in zip(mi[keep], mj[keep]):
            adj.setdefault(int(i), set()).add(int(j))
            adj.setdefault(int(j), set()).add(int(i))
    return adj


def grow_connected(adj, seed, n, rng):
    """Grow a connected n-molecule subgraph from ``seed``.

    Random frontier expansion: at each step a molecule adjacent to the current
    cluster is chosen uniformly and attached through a random one of its
    neighbours already inside. Uniform choice from the frontier is a cheap
    approximation to uniform connected-subgraph sampling; it mildly favours
    high-degree regions.

    Returns
    -------
    (order, parent) or None
        ``order`` is the list of molecule indices in growth order; ``parent``
        maps each to the molecule it attached through (None for the seed).
        None if the connected component is smaller than n.
    """
    order = [seed]
    parent = {seed: None}
    inside = {seed}
    frontier = set(adj.get(seed, ()))
    while len(order) < n:
        if not frontier:
            return None
        cands = sorted(frontier)
        c = int(cands[rng.integers(len(cands))])
        anchors = sorted(adj[c] & inside)
        parent[c] = int(anchors[rng.integers(len(anchors))])
        order.append(c)
        inside.add(c)
        frontier |= adj.get(c, set())
        frontier -= inside
    return order, parent


def unwrap_shifts(order, parent, frac_centroids, periodic):
    """Integer lattice shifts making each molecule contiguous with its parent.

    Each molecule is shifted to the image nearest the neighbour it attached
    through, so the whole cluster is contiguous in space.
    """
    shifts = {}
    for m in order:
        if parent[m] is None or not periodic:
            shifts[m] = np.zeros(3)
        else:
            p = parent[m]
            u_p = frac_centroids[p] + shifts[p]
            shifts[m] = np.round(u_p - frac_centroids[m])
    return shifts


def _spreads(cen):
    """Radius of gyration and largest pairwise distance of a set of points."""
    cen0 = cen.mean(axis=0)
    rg = float(np.sqrt(((cen - cen0) ** 2).sum(axis=1).mean()))
    n = len(cen)
    dmax = float(
        max(np.linalg.norm(cen[i] - cen[j]) for i in range(n) for j in range(i))
    )
    return rg, dmax


def pilot_spreads(
    adj, n, frac_centroids, T, periodic, spread_metric, n_mol, rng, n_pilot
):
    """Spread values of an unstratified pilot sample, used to place bin edges."""
    out = []
    tries = 0
    while len(out) < n_pilot and tries < 20 * n_pilot:
        tries += 1
        g = grow_connected(adj, int(rng.integers(n_mol)), n, rng)
        if g is None:
            continue
        order, parent = g
        shifts = unwrap_shifts(order, parent, frac_centroids, periodic)
        cen = np.array([(frac_centroids[m] + shifts[m]) @ T for m in order])
        rg, dmax = _spreads(cen)
        out.append(rg if spread_metric == "rg" else dmax)
    return np.array(out)


def bond_index_pairs(configuration):
    """The bonds of a configuration as 0-based atom-index pairs, with orders.

    For a configuration without symmetry operators (the normal case, and every
    non-periodic one) the pairs come straight from the bond table, since
    ``symmetry.bond_atoms`` is only populated for periodic configurations. For
    a symmetric crystal the symmetry-expanded list is used.

    Returns
    -------
    (pairs, orders)
        ``pairs`` is a list of (i, j) index tuples; ``orders`` the matching
        list of bond orders, or None if it could not be matched up.
    """
    atoms = configuration.atoms
    bonds = configuration.bonds
    if configuration.symmetry.n_symops == 1:
        index = {aid: k for k, aid in enumerate(atoms.ids)}
        Is = bonds.get_column_data("i")
        Js = bonds.get_column_data("j")
        pairs = [(index[i], index[j]) for i, j in zip(Is, Js)]
        orders = list(bonds.get_column_data("bondorder"))
    else:
        pairs = [(int(i), int(j)) for i, j in configuration.symmetry.bond_atoms]
        orders = list(bonds.bondorders)
    if len(orders) != len(pairs):
        orders = None
    return pairs, orders


def make_put_property(conf, tag=PROPERTY_TAG):
    """Return a closure storing ``name+tag`` properties on ``conf``.

    Defines the property first if the database lacks it (the properties are
    normally pre-registered from ``data/properties.csv``).
    """
    props = conf.properties

    def _put(name, value, units=None, _type="float"):
        full = f"{name}{tag}"
        if not props.exists(full):
            props.add(full, _type, units=units, noerror=True)
        props.put(full, value)

    return _put


def extract_nmers(
    configuration,
    n,
    n_samples,
    *,
    cutoff=3.5,
    contact_elements=None,
    spread_metric="rg",
    spread_bins=None,
    balance_motifs=False,
    system=None,
    system_name="clusters",
    name_prefix="",
    max_attempts=None,
    rng=None,
    store_properties=True,
    seen=None,
):
    """Extract n-molecule clusters from a (periodic) configuration.

    The clusters are returned as unwrapped, non-periodic configurations
    suitable for cluster QM.

    Selection is by *connected subgraph* of a molecular contact graph, not by
    nearest neighbours: every molecule in a returned n-mer is within
    ``cutoff`` of at least one other, but chains, rings and stars all occur,
    and seeds are random so a single frame yields many distinct clusters.

    Parameters
    ----------
    configuration : molsystem _Configuration
        Source configuration (periodic or not). Molecules are found by bonds.
    n : int
        Molecules per cluster (>= 2).
    n_samples : int
        Target number of clusters (fewer if the pool or quotas run out).
    cutoff : float
        Contact distance (Å) defining the neighbour graph. 3.5 with
        ``contact_elements=["O"]`` is the usual H-bond criterion for water.
    contact_elements : list[str] or None
        Elements whose atoms define contact; None uses all atoms.
    spread_metric : "rg" | "dmax"
        Compactness coordinate used for stratification: radius of gyration of
        the molecular centroids, or the largest centroid-centroid distance.
    spread_bins : int, sequence of float, or None
        Bin edges on the spread coordinate. An int asks for that many
        equal-quantile bins, with edges taken from an unstratified pilot
        sample of the same frame -- the usual choice, since sensible edges
        depend on n and on the system. If given, acceptance is quota-
        limited so bins fill roughly evenly (flat-in-spread, analogous to the
        dimer builder's flat-in-energy selection). Samples outside the edges
        are rejected. None disables stratification.
    balance_motifs : bool
        Also balance across contact-graph topology (chain/ring/star/...),
        jointly with the spread bins.
    system : molsystem _System or None
        Destination system. If None, ``system_name`` is looked up in the
        configuration's database and created if needed.
    system_name : str
        Destination system name, used when ``system`` is None.
    name_prefix : str
        Prefix (e.g. a frame label) for the configuration names, which are
        ``f"{name_prefix}{seed}_{m1-m2-...}"`` and unique within a frame.
    max_attempts : int or None
        Sampling attempts before giving up (default 50 * n_samples).
    rng : int, numpy Generator or None
        Seed / generator for reproducibility.
    store_properties : bool
        Store size, spread, motif, edge count, and source molecule list as
        ``*#ExtractClusters#scan`` properties on each configuration.
    seen : set or None
        Molecule sets already emitted (as frozensets); updated in place. Pass
        the same set across calls on the same frame to avoid duplicates.

    Returns
    -------
    (configurations, records, info)
        New configurations; one metadata dict per cluster with keys
        ``name, n, molecules, seed, order, rg, dmax, n_edges, degrees, motif,
        bin``; and an ``info`` dict with ``edges`` (the bin edges used, or
        None), ``attempts``, ``per_key`` quota, ``n_molecules`` and
        ``n_contacts`` of the frame's graph, and ``warning`` (a message or
        None) if fewer than requested were produced.
    """
    if n < 2:
        raise ValueError("The cluster size must be >= 2")
    if spread_metric not in SPREAD_METRICS:
        raise ValueError(f"spread_metric must be one of {SPREAD_METRICS}")
    if n_samples < 1:
        raise ValueError("The number of clusters must be >= 1")
    rng = np.random.default_rng(rng)
    if max_attempts is None:
        max_attempts = 50 * n_samples
    if seen is None:
        seen = set()

    periodic = configuration.periodicity != 0
    molecules = [np.asarray(m) for m in configuration.find_molecules(as_indices=True)]
    n_mol = len(molecules)
    if n_mol < n:
        raise ValueError(
            f"The configuration has only {n_mol} molecules; cannot make {n}-mers"
        )
    atom_mol = np.empty(configuration.n_atoms, dtype=int)
    for k, m in enumerate(molecules):
        atom_mol[m] = k

    if periodic:
        # Fractional, molecules kept whole across the boundary.
        frac = configuration.atoms.get_coordinates(
            fractionals=True, in_cell="molecule", as_array=True
        )
        T = configuration.cell.to_cartesians_transform(as_array=True)
        cell_parameters = configuration.cell.parameters
    else:
        frac = configuration.atoms.get_coordinates(fractionals=False, as_array=True)
        T = np.eye(3)
        cell_parameters = None

    frac_centroids = np.array([frac[m].mean(axis=0) for m in molecules])

    symbols = np.array(configuration.atoms.symbols)
    if contact_elements is None:
        contact_mask = np.ones(len(symbols), dtype=bool)
    else:
        contact_mask = np.isin(symbols, list(contact_elements))
        if not contact_mask.any():
            raise ValueError(
                f"No atoms match the contact elements {list(contact_elements)}"
            )

    adj = molecule_neighbour_graph(
        frac, atom_mol, contact_mask, T, cutoff, periodic, cell_parameters
    )
    if not adj:
        raise ValueError(f"No molecular contacts within {cutoff} Å")
    n_contacts = sum(len(v) for v in adj.values()) // 2

    # Stratification bins -----------------------------------------------------
    if isinstance(spread_bins, (int, np.integer)):
        # Choose edges from the data: equal-quantile bins of the spread
        # coordinate over an unstratified pilot sample.
        n_bins_requested = int(spread_bins)
        pilot = pilot_spreads(
            adj,
            n,
            frac_centroids,
            T,
            periodic,
            spread_metric,
            n_mol,
            rng,
            n_pilot=max(200, 20 * n_bins_requested),
        )
        if len(pilot) < 2 * n_bins_requested:
            raise ValueError(
                f"Too few {n}-mers found ({len(pilot)}) to choose "
                f"{n_bins_requested} spread bins"
            )
        edges = np.quantile(pilot, np.linspace(0.0, 1.0, n_bins_requested + 1))
        edges[0] -= 1e-9
        edges[-1] += 1e-9
    elif spread_bins is None:
        edges = None
    else:
        edges = np.asarray(spread_bins, dtype=float)
        if edges.ndim != 1 or len(edges) < 2:
            raise ValueError("spread_bins must give at least two bin edges")
        if np.any(np.diff(edges) <= 0):
            raise ValueError("spread_bins edges must be strictly increasing")
    n_bins = 1 if edges is None else len(edges) - 1
    stratified = edges is not None or balance_motifs

    atnos = np.asarray(configuration.atoms.atomic_numbers)
    bond_pairs, bond_orders = bond_index_pairs(configuration)
    bond_indices = {}
    for k, (i, j) in enumerate(bond_pairs):
        bond_indices.setdefault(int(i), []).append(k)

    if system is None:
        system_db = configuration.system_db
        systems = system_db.get_systems(system_name)
        system = systems[0] if systems else system_db.create_system(name=system_name)

    # Quotas: each (motif?, bin) key gets an equal share of n_samples over the
    # keys seen so far -- the motif set is not known in advance.
    per_key = None
    configurations = []
    records = []
    counts = Counter()
    attempts = 0
    while len(configurations) < n_samples and attempts < max_attempts:
        attempts += 1
        seed = int(rng.integers(n_mol))
        grown = grow_connected(adj, seed, n, rng)
        if grown is None:
            continue
        order, parent = grown
        key = frozenset(order)
        if key in seen:
            continue

        shifts = unwrap_shifts(order, parent, frac_centroids, periodic)
        cen = np.array([(frac_centroids[m] + shifts[m]) @ T for m in order])
        rg, dmax = _spreads(cen)
        spread = rg if spread_metric == "rg" else dmax

        inside = set(order)
        sub_edges = [
            (i, j) for i in order for j in adj.get(i, ()) if j in inside and i < j
        ]
        n_edges = len(sub_edges)
        deg = Counter()
        for i, j in sub_edges:
            deg[i] += 1
            deg[j] += 1
        degrees = sorted((deg[m] for m in order), reverse=True)
        motif = classify_motif(n, n_edges, degrees)

        if edges is not None:
            b = int(np.digitize(spread, edges)) - 1
            if b < 0 or b >= n_bins:
                continue  # outside requested range
        else:
            b = 0
        qkey = (motif, b) if balance_motifs else (b,)
        if stratified:
            # equal share per key, over bins x (motifs seen so far)
            if balance_motifs:
                n_motifs = max(1, len({k[0] for k in counts} | {motif}))
            else:
                n_motifs = 1
            per_key = math.ceil(n_samples / (n_bins * n_motifs))
            if counts[qkey] >= per_key:
                continue

        seen.add(key)
        counts[qkey] += 1

        # Build the cluster ---------------------------------------------------
        atom_indices = np.concatenate([molecules[m] for m in order])
        xyz = np.concatenate([(frac[molecules[m]] + shifts[m]) @ T for m in order])
        xyz -= xyz.mean(axis=0)  # centre the cluster at the origin
        to_new = {int(old): new for new, old in enumerate(atom_indices)}
        Is, Js, orders = [], [], []
        for i in atom_indices:
            for k in bond_indices.get(int(i), ()):
                i0, j0 = bond_pairs[k]
                if j0 in to_new:
                    Is.append(to_new[int(i0)])
                    Js.append(to_new[int(j0)])
                    if bond_orders is not None:
                        orders.append(bond_orders[k])

        name = f"{name_prefix}{seed}_" + "-".join(str(m) for m in sorted(order))
        conf = system.create_configuration(
            name=name, periodicity=0, coordinate_system="Cartesian", make_current=False
        )
        ids = conf.atoms.append(
            atno=atnos[atom_indices].tolist(),
            x=xyz[:, 0].tolist(),
            y=xyz[:, 1].tolist(),
            z=xyz[:, 2].tolist(),
        )
        if Is:
            kwargs = {}
            if bond_orders is not None:
                kwargs["bondorder"] = orders
            conf.bonds.append(i=[ids[k] for k in Is], j=[ids[k] for k in Js], **kwargs)

        rec = dict(
            name=name,
            n=n,
            molecules=sorted(order),
            seed=seed,
            order=list(order),
            rg=rg,
            dmax=dmax,
            n_edges=n_edges,
            degrees=degrees,
            motif=motif,
            bin=b,
        )
        if store_properties:
            _put = make_put_property(conf)
            _put("cluster size", n, None, "int")
            _put("spread rg", rg, "Å")
            _put("spread dmax", dmax, "Å")
            _put("n edges", n_edges, None, "int")
            _put("motif", motif, None, "str")
            _put("spread bin", b, None, "int")
            _put(
                "source molecules", "-".join(str(m) for m in sorted(order)), None, "str"
            )
        configurations.append(conf)
        records.append(rec)

    warning = None
    if len(configurations) < n_samples:
        warning = (
            f"Produced {len(configurations)} of {n_samples} requested {n}-mers "
            f"after {attempts} attempts"
        )
        if stratified:
            warning += (
                f" (per-key quota {per_key}). Quotas for (motif, bin) cells that "
                "never occur cannot fill -- e.g. rings only exist in the compact "
                "bins -- so a short fill here usually reflects the physics rather "
                "than a sampling problem. Check the motif x bin table and revisit "
                "the stratification or motif balancing if needed."
            )
        else:
            warning += (
                ". The frame may not contain that many distinct connected "
                f"{n}-mers within the contact cutoff."
            )

    info = dict(
        edges=None if edges is None else [float(e) for e in edges],
        attempts=attempts,
        per_key=per_key,
        n_molecules=n_mol,
        n_contacts=n_contacts,
        warning=warning,
    )
    return configurations, records, info


def cluster_summary(records, edges=None):
    """Tabulate how a set of extracted n-mers is distributed over motif and
    spread bin -- the quick check that the set is actually balanced.

    Parameters
    ----------
    records : [dict]
        The records from :func:`extract_nmers`.
    edges : [float] or None
        The bin edges, printed as a final line if given.

    Returns
    -------
    str
    """
    tab = Counter((r["motif"], r["bin"]) for r in records)
    motifs = sorted({m for m, _ in tab})
    bins = sorted({b for _, b in tab})
    w = max(8, max((len(m) for m in motifs), default=8))
    header = " " * w + "".join(f"{('bin%d' % b):>8}" for b in bins) + f"{'total':>8}"
    lines = [header]
    for m in motifs:
        row = [tab[(m, b)] for b in bins]
        lines.append(f"{m:<{w}}" + "".join(f"{c:>8}" for c in row) + f"{sum(row):>8}")
    tot = [sum(tab[(m, b)] for m in motifs) for b in bins]
    lines.append(f"{'total':<{w}}" + "".join(f"{c:>8}" for c in tot) + f"{sum(tot):>8}")
    if edges is not None:
        e = list(edges)
        lines.append(
            "bins (Å): "
            + ", ".join(f"{e[i]:.2f}-{e[i + 1]:.2f}" for i in range(len(e) - 1))
        )
    return "\n".join(lines)
