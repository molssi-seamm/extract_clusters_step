===============
Developer Guide
===============

.. toctree::
   :maxdepth: 1

   campaigns/2026-09-17/index
   campaigns/2026-09-18/index
   modules
   ../installation
   ../contributing
   ../usage

Code layout
-----------

``extract_clusters_step/cluster_sampling.py``
    The algorithm, free of any SEAMM dependency: ``molecule_neighbour_graph``,
    ``grow_connected``, ``unwrap_shifts``, ``classify_motif``,
    ``bond_index_pairs``, ``extract_nmers`` and ``cluster_summary``. Unit
    tested directly in ``tests/test_cluster_sampling.py`` (pure numpy) and
    ``tests/test_extract_nmers.py`` (on a synthetic periodic water box).
``extract_clusters_step/extract_clusters.py``
    The ``seamm.Node``: parses the parameters, calls ``extract_nmers`` once
    per requested size sharing one ``seen`` set, writes ``clusters.csv`` and
    prints the motif × bin summary.
``extract_clusters_step/tk_extract_clusters.py``
    The dialog; ``reset_dialog`` shows only the stratification controls that
    apply to the chosen scheme.
``extract_clusters_step/data/properties.csv``
    Registers the ``#ExtractClusters#scan`` properties with molsystem at
    import.
