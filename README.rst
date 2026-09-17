==============================
SEAMM Extract Clusters Plug-in
==============================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/extract_clusters_step
   :target: https://github.com/molssi-seamm/extract_clusters_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/extract_clusters_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/extract_clusters_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/extract_clusters_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/extract_clusters_step
   :alt: Code Coverage

.. image:: https://img.shields.io/lgtm/grade/python/g/molssi-seamm/extract_clusters_step.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/molssi-seamm/extract_clusters_step/context:python
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/extract_clusters_step/workflows/Documentation/badge.svg
   :target: https://molssi-seamm.github.io/extract_clusters_step/index.html
   :alt: Documentation Status

.. image:: https://pyup.io/repos/github/molssi-seamm/extract_clusters_step/shield.svg
   :target: https://pyup.io/repos/github/molssi-seamm/extract_clusters_step/
   :alt: Updates for Dependencies

.. image:: https://img.shields.io/pypi/v/extract_clusters_step.svg
   :target: https://pypi.python.org/pypi/extract_clusters_step
   :alt: PyPi VERSION

A SEAMM plug-in for extracting molecular clusters (trimers, tetramers, ... larger
n-mers) from a condensed-phase, typically periodic, configuration as unwrapped,
non-periodic structures -- e.g. many-body training data and diagnostics for
machine-learned force fields.

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/extract_clusters_step/index.html
* Code: https://github.com/molssi-seamm/extract_clusters_step

Features
--------

* Clusters are **connected subgraphs of a molecular contact graph** (molecules
  are in contact if their contact atoms are within a cutoff, minimum image), so
  chains, rings and stars all occur -- not just the most compact cluster.
* Any cluster size, and several sizes per frame (e.g. ``3, 4``).
* **Stratification** so the set is flat in a spread coordinate (radius of
  gyration or largest centroid separation), with bin edges from equal quantiles
  of a pilot sample or given explicitly; optional balancing over the
  contact-graph motif.
* Clusters are unwrapped across the periodic boundary, centred, non-periodic,
  with molecules and bonds intact.
* Provenance on every configuration: unique names
  ``<frame>_<seed>_<m1-m2-...>`` and ``#ExtractClusters#scan`` properties (size,
  spread, motif, contacts, source molecules) that survive SDF/extxyz export;
  a ``clusters.csv`` summary per step.
* Works on any cell (orthorhombic fast path via a periodic KD-tree; exact
  minimum image otherwise) and on non-periodic sources.

Acknowledgements
----------------

This package was created with Cookiecutter_ and the
`molssi-seamm/cookiecutter-seamm-plugin`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award ACI-1547580

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
