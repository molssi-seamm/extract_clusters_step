=====================================================
2026-09-17: Cluster (n-mer) extraction for MLFF data
=====================================================

This is the design record for the plug-in, moved here from the workspace
planning document (``~/Work/SEAMM/nmer-extractor-plan.md``) once the package
existed. The science context is canonical in the ``~/Sites`` lab notebook:

- ``~/Sites/mlff-training/2026-09-15_water-liquid-density-diagnostics/``
- ``~/Sites/training-data/2026-08-13_h2o-trimer-run1/``
- ``~/Sites/mlff-training/2026-09-16_trimer-mlff-liquid-density/``

Why
---

A water MLFF trained on dimers only ran NPT at 0.664 g/cm³ with hydrogen bonds
frozen at the gas-phase dimer distance. Adding a trimer set recovered the
liquid (1.037 g/cm³, correct O–O peak and coordination); at the experimental
density the remaining defect is the equation of state alone, ≈1.7 kJ/mol per
molecule of excess cohesion -- the size of the omitted four-body term. Hence a
general way to pull n-mers (n = 3, 4, 5, and larger clusters) out of
condensed-phase configurations, to measure E₄ directly and to add balanced
n-mer sets to the training data.

Design decisions
----------------

1. **Connected-subgraph sampling, not nearest neighbours.** Nearest-neighbour
   selection gives one deterministic, maximally compact cluster per frame and
   kills topological variety; for water the linear hydrogen-bond chain is
   exactly where cooperativity is largest. An n-mer is a connected subgraph of
   the molecular contact graph (A~B if any contact atom of A is within the
   cutoff of one of B, minimum image), grown by random frontier expansion from
   a random seed.
2. **Stratify, don't just accept.** Flat-in-spread acceptance (radius of
   gyration or largest centroid separation) with quantile edges from a pilot
   sample -- the dimer builder's flat-in-energy philosophy -- plus optional
   joint balancing on motif.
3. **Gas-phase clusters done right:** unwrap by parent (each molecule made
   contiguous with the neighbour it attached through), centre, periodicity 0,
   molecules and bonds intact.
4. **Provenance on every configuration:** unique names
   ``{prefix}{seed}_{m1-m2-...}`` (SEAMM loops key on configuration names) and
   ``*#ExtractClusters#scan`` properties that survive SDF export.
5. **General, not water-specific:** molecules by bonds, contact criterion by
   element list or all atoms, any cell (periodic KD-tree for orthorhombic,
   exact chunked minimum image otherwise), non-periodic sources too.

Verified
--------

- Both contact-graph paths agree with brute-force minimum image; grown
  clusters are connected and contiguous after unwrapping; motif labels are
  correct (``tests/test_cluster_sampling.py``).
- End to end on a synthetic periodic box and on the real 1667-water frame
  (``~/Work/SEAMM/Testing/extract_clusters.flow`` + ``water_box.cif``): 30
  trimers + 30 tetramers, flat over 3 quantile bins, unique names, 3n atoms /
  2n bonds, centred, connected, properties present in the written SDF.

Known behaviour
---------------

- **Topology and spread are correlated.** With motif balancing, (ring, outer
  bin) quota cells can never fill and the call returns short with a warning
  naming the cause. The quota policy (redistribute, balance one axis, or accept
  unequal cells) is an open design choice; do not "fix" it silently -- it is
  real physics.
- Hand-picked bin edges are usually wrong (no water trimer has rg < 1.6 Å);
  prefer quantile bins.
- molsystem: ``from_mmcif_text`` drops the ``_cell`` tags of a SEAMM-written
  mmCIF and loads it non-periodic (use ``.cif`` for periodic frames);
  ``symmetry.bond_atoms`` is never populated for periodicity 0, so
  ``bond_index_pairs`` reads the bond table directly when there is no
  symmetry.

Open items
----------

1. Quota redistribution vs one-axis balancing (above).
2. Seed-environment stratification (coordination number / tetrahedral order
   q) so minority environments are sampled.
3. Large clusters (n ≈ 10–30) to fill the MLFF receptive field: same sampler,
   motif naming beyond n = 4 is ``e<edges>``.
4. Mixed systems (EC/FEC + Li⁺/BF₄⁻): per-element-pair contact criteria;
   charged fragments need net charge/multiplicity carried to the CP step.
5. Performance for ≫10⁴ atoms with non-orthorhombic cells (general path is
   O(N²) on contact atoms; a supercell KD-tree would fix it).
6. A trajectory-frame loader into a periodic molsystem configuration.
7. First real use: tetramers from the ρ = 0.997 NVT trajectory → N-fragment
   counterpoise labelling → MLFF vs revDSD tetramer energies (direct E₄
   diagnostic).
