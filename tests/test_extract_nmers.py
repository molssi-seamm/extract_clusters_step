# -*- coding: utf-8 -*-

"""End-to-end tests of extract_nmers on a periodic water box."""

import numpy as np
import pytest

import extract_clusters_step
from extract_clusters_step.cluster_sampling import (
    PROPERTY_TAG,
    check_bonds,
    cluster_summary,
    extract_nmers,
)


def _prop(conf, name):
    full = f"{name}{PROPERTY_TAG}"
    return conf.properties.get(full, match="exact")[full]["value"]


def _check_clusters(confs, records, n, cutoff=3.5, contact="O"):
    assert len(confs) == len(records)
    for c, r in zip(confs, records):
        assert c.periodicity == 0
        assert c.n_atoms == 3 * n
        assert c.bonds.n_bonds == 2 * n
        assert set(c.bonds.bondorders) == {1}
        xyz = c.atoms.get_coordinates(fractionals=False, as_array=True)
        sym = np.array(c.atoms.symbols)
        # centred at the origin
        assert np.abs(xyz.mean(0)).max() < 1e-6
        # connected under plain (unwrapped) distances
        if contact == "all":
            pts = xyz
            mol = np.repeat(np.arange(n), 3)
        else:
            pts = xyz[sym == contact]
            mol = np.arange(n)
        close = np.linalg.norm(pts[:, None] - pts[None], axis=2) < cutoff
        A = np.zeros((n, n), dtype=bool)
        for a, b in zip(*np.nonzero(close)):
            if mol[a] != mol[b]:
                A[mol[a], mol[b]] = True
        seen = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in np.nonzero(A[u])[0]:
                if v not in seen:
                    seen.add(int(v))
                    stack.append(int(v))
        assert len(seen) == n
        # molecules intact
        for k in range(n):
            o, h1, h2 = xyz[3 * k], xyz[3 * k + 1], xyz[3 * k + 2]
            assert 0.8 < np.linalg.norm(h1 - o) < 1.2
            assert 0.8 < np.linalg.norm(h2 - o) < 1.2
        # stored properties match the record
        assert abs(_prop(c, "spread rg") - r["rg"]) < 1e-9
        assert abs(_prop(c, "spread dmax") - r["dmax"]) < 1e-9
        assert _prop(c, "motif") == r["motif"]
        assert _prop(c, "cluster size") == n
        assert _prop(c, "n edges") == r["n_edges"]
        assert _prop(c, "spread bin") == r["bin"]
        assert _prop(c, "source molecules") == "-".join(str(m) for m in r["molecules"])
        assert c.name == r["name"]


def test_box_fixture(water_box):
    assert water_box.periodicity == 3
    assert water_box.n_atoms == 3 * 216
    assert len(water_box.find_molecules()) == 216


def test_trimers_stratified(water_box):
    confs, recs, info = extract_nmers(
        water_box,
        3,
        30,
        cutoff=3.5,
        contact_elements=["O"],
        spread_bins=3,
        system_name="trimers",
        name_prefix="f0_",
        rng=7,
    )
    assert len(confs) == 30
    assert info["warning"] is None
    assert len(info["edges"]) == 4
    _check_clusters(confs, recs, 3)
    names = [c.name for c in confs]
    assert len(names) == len(set(names))
    assert all(name.startswith("f0_") for name in names)
    # flat in spread: every bin has its equal share
    bins = [r["bin"] for r in recs]
    assert sorted(set(bins)) == [0, 1, 2]
    assert max(bins.count(b) for b in (0, 1, 2)) == 10
    # motifs are the trimer ones
    assert set(r["motif"] for r in recs) <= {"chain", "ring"}
    # summary table has a row per motif and a total line
    table = cluster_summary(recs, info["edges"])
    assert "total" in table and "bins (Å)" in table
    # the destination system was created
    assert water_box.system_db.get_systems("trimers")[0].n_configurations == 30


def test_tetramers_unstratified_and_dedupe(water_box):
    seen = set()
    c1, r1, i1 = extract_nmers(
        water_box, 4, 20, contact_elements=["O"], system_name="tets", rng=1, seen=seen
    )
    c2, r2, i2 = extract_nmers(
        water_box, 4, 20, contact_elements=["O"], system_name="tets", rng=2, seen=seen
    )
    _check_clusters(c1 + c2, r1 + r2, 4)
    assert i1["edges"] is None and i1["per_key"] is None
    keys = [frozenset(r["molecules"]) for r in r1 + r2]
    assert len(keys) == len(set(keys)) == 40
    assert set(r["motif"] for r in r1 + r2) <= {
        "chain",
        "star",
        "ring",
        "paw",
        "diamond",
        "K4",
    }


def test_reproducible_with_seed(water_box):
    _, r1, _ = extract_nmers(water_box, 3, 10, contact_elements=["O"], rng=3)
    _, r2, _ = extract_nmers(water_box, 3, 10, contact_elements=["O"], rng=3)
    assert [r["molecules"] for r in r1] == [r["molecules"] for r in r2]


def test_explicit_edges_and_short_fill(water_box):
    # Edges that exclude everything except a narrow band -> short fill, warning
    confs, recs, info = extract_nmers(
        water_box,
        3,
        50,
        contact_elements=["O"],
        spread_bins=[1.0, 1.05],
        max_attempts=200,
        rng=0,
    )
    assert len(confs) < 50
    assert info["warning"] is not None
    assert "quota" in info["warning"]
    for r in recs:
        assert 1.0 <= r["rg"] < 1.05


def test_balance_motifs(water_box):
    confs, recs, info = extract_nmers(
        water_box, 3, 40, contact_elements=["O"], balance_motifs=True, rng=11
    )
    motifs = {r["motif"] for r in recs}
    # both trimer motifs are present in a dense box and each is capped at its share
    if len(motifs) == 2:
        counts = [sum(1 for r in recs if r["motif"] == m) for m in motifs]
        assert max(counts) <= 20


def test_dmax_metric_and_all_atoms(water_box):
    confs, recs, info = extract_nmers(
        water_box, 3, 12, cutoff=2.5, spread_metric="dmax", spread_bins=2, rng=5
    )
    assert len(confs) == 12
    # With all-atom contacts at 2.5 Å the O-O distances can exceed 3.5 Å, so
    # check connectivity with the same all-atom criterion.
    _check_clusters(confs, recs, 3, cutoff=2.5, contact="all")


def test_non_periodic_source(water_box):
    """A non-periodic source (a cluster) works too; no unwrapping is needed."""
    confs, recs, info = extract_nmers(
        water_box, 6, 1, contact_elements=["O"], system_name="hexamers", rng=4
    )
    assert len(confs) == 1
    hexamer = confs[0]
    sub, subr, _ = extract_nmers(
        hexamer, 3, 5, contact_elements=["O"], system_name="sub", rng=4
    )
    assert 1 <= len(sub) <= 5
    _check_clusters(sub, subr, 3)


def test_errors(water_box):
    with pytest.raises(ValueError):
        extract_nmers(water_box, 1, 5)
    with pytest.raises(ValueError):
        extract_nmers(water_box, 3, 5, spread_metric="bogus")
    with pytest.raises(ValueError):
        extract_nmers(water_box, 3, 5, contact_elements=["Xx"])
    with pytest.raises(ValueError):
        extract_nmers(water_box, 3, 5, cutoff=0.1, contact_elements=["O"])
    with pytest.raises(ValueError):
        extract_nmers(water_box, 3, 5, spread_bins=[2.0, 1.0])
    with pytest.raises(ValueError):
        extract_nmers(water_box, 1000, 5)


def test_sdf_round_trip(water_box):
    confs, recs, _ = extract_nmers(
        water_box, 3, 2, contact_elements=["O"], name_prefix="f0_", rng=9
    )
    txt = confs[0].to_sdf_text()
    assert txt.splitlines()[0].strip().endswith(confs[0].name)
    assert f"motif{PROPERTY_TAG}" in txt
    assert f"source molecules{PROPERTY_TAG}" in txt


def test_node_run(water_box, tmp_path, monkeypatch):
    """Run the node itself against the water box, without a flowchart."""
    import seamm

    node = extract_clusters_step.ExtractClusters()
    node._id = (1,)
    monkeypatch.setattr(type(node), "directory", property(lambda self: str(tmp_path)))
    db = water_box.system_db
    monkeypatch.setattr(node, "get_variable", lambda name: db)
    monkeypatch.setattr(
        node, "get_system_configuration", lambda P=None: (water_box.system, water_box)
    )
    monkeypatch.setattr(seamm.Node, "run", lambda self, printer=None: None)
    if seamm.flowchart_variables is None:
        monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())

    P = node.parameters
    P["cluster sizes"].value = "3, 4"
    P["number of clusters"].value = 12
    P["contact elements"].value = "O"
    P["number of bins"].value = 2
    P["balance motifs"].value = "yes"
    P["random seed"].value = "3"
    P["make current"].value = "no"
    assert node.run() is None

    clusters = db.get_systems("water box clusters")[0]
    assert clusters.n_configurations > 0
    names = [c.name for c in clusters.configurations]
    assert all(name.startswith("frame0_") for name in names)
    assert len(names) == len(set(names))
    # properties were stored (the boolean came through as True at run time)
    c0 = clusters.configurations[0]
    assert _prop(c0, "cluster size") in (3, 4)
    # 'make current: no' left the source system current
    assert db.system.name == "water box"
    assert (tmp_path / "clusters.csv").exists()
    assert len((tmp_path / "clusters.csv").read_text().splitlines()) == len(names) + 1

    # and 'make current: yes' switches to the clusters
    P["make current"].value = "yes"
    P["random seed"].value = "4"
    node.run()
    assert db.system.name == "water box clusters"


def test_bondless_configuration_is_rejected(water_box):
    """No bonds + molecular elements -> a clear error, not atom 'clusters'."""
    water_box.bonds.delete()
    assert water_box.bonds.n_bonds == 0
    with pytest.raises(ValueError, match="Perceive bonds"):
        extract_nmers(water_box, 3, 5, contact_elements=["O"])
    # perceiving the bonds makes it work again (molsystem >= 2026.9.17)
    if not hasattr(water_box, "perceive_bonds"):
        pytest.skip("molsystem without perceive_bonds")
    assert water_box.perceive_bonds() == 2 * 216
    confs, recs, info = extract_nmers(water_box, 3, 5, contact_elements=["O"], rng=1)
    assert len(confs) == 5


def test_bondless_rare_gas_is_allowed(db):
    """A bond-less argon box is legitimately one atom per molecule."""
    rng = np.random.default_rng(0)
    L = 20.0
    xyz = rng.uniform(0, L, (100, 3))
    conf = db.create_system(name="Ar").create_configuration(
        name="liquid",
        periodicity=3,
        coordinate_system="Cartesian",
        cell_parameters=[L, L, L, 90, 90, 90],
    )
    conf.atoms.append(
        atno=[18] * 100,
        x=xyz[:, 0].tolist(),
        y=xyz[:, 1].tolist(),
        z=xyz[:, 2].tolist(),
    )
    check_bonds(conf)  # no error
    confs, recs, info = extract_nmers(conf, 3, 5, cutoff=4.5, rng=2)
    assert len(confs) == 5
    assert all(c.n_atoms == 3 for c in confs)
