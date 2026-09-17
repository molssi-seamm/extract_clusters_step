# -*- coding: utf-8 -*-

"""Fixtures: a small periodic water box as a molsystem configuration."""

import math

import numpy as np
import pytest

from molsystem.system_db import SystemDB

_db_counter = [0]


def water_box_configuration(db, n_side=6, spacing=3.1, name="frame0"):
    """A periodic box of waters on a jittered cubic lattice.

    ``spacing`` ~3.1 Å puts the O-O nearest neighbours inside the 3.5 Å
    H-bond cutoff so the contact graph is well connected.
    """
    rng = np.random.default_rng(12345)
    L = n_side * spacing
    r0 = 0.9572
    theta0 = math.radians(104.52)
    xh = r0 * math.sin(theta0 / 2)
    zh = r0 * math.cos(theta0 / 2)
    xs, ys, zs, atnos = [], [], [], []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                o = np.array([i, j, k]) * spacing + rng.uniform(-0.4, 0.4, 3)
                # random orientation
                q = rng.normal(size=4)
                q /= np.linalg.norm(q)
                a, b, c, d = q
                R = np.array(
                    [
                        [
                            a * a + b * b - c * c - d * d,
                            2 * (b * c - a * d),
                            2 * (b * d + a * c),
                        ],
                        [
                            2 * (b * c + a * d),
                            a * a - b * b + c * c - d * d,
                            2 * (c * d - a * b),
                        ],
                        [
                            2 * (b * d - a * c),
                            2 * (c * d + a * b),
                            a * a - b * b - c * c + d * d,
                        ],
                    ]
                )
                for v in ([0, 0, 0], [xh, 0, zh], [-xh, 0, zh]):
                    p = o + R @ np.array(v)
                    xs.append(float(p[0]))
                    ys.append(float(p[1]))
                    zs.append(float(p[2]))
                atnos.extend([8, 1, 1])

    system = db.create_system(name="water box")
    conf = system.create_configuration(
        name=name,
        periodicity=3,
        coordinate_system="Cartesian",
        cell_parameters=[L, L, L, 90, 90, 90],
        make_current=True,
    )
    ids = conf.atoms.append(atno=atnos, x=xs, y=ys, z=zs)
    n_mol = n_side**3
    Is = [ids[3 * m] for m in range(n_mol) for _ in (0, 1)]
    Js = [ids[3 * m + h] for m in range(n_mol) for h in (1, 2)]
    conf.bonds.append(i=Is, j=Js, bondorder=[1] * len(Is))
    return conf


@pytest.fixture()
def db():
    _db_counter[0] += 1
    db = SystemDB(
        filename=f"file:extract_test_{_db_counter[0]}?mode=memory&cache=shared"
    )
    yield db
    db.close()


@pytest.fixture()
def water_box(db):
    return water_box_configuration(db)
