=======
History
=======
2026.9.17 -- Initial release of the Extract Clusters step
    * Extracts n-molecule clusters (trimers, tetramers, ... larger n-mers) from the
      current, typically periodic, condensed-phase configuration as unwrapped,
      centred, non-periodic configurations in a new system, e.g. many-body training
      data and diagnostics for machine-learned force fields.
    * Clusters are connected subgraphs of a molecular contact graph (molecules are in
      contact if their contact atoms are within a cutoff, minimum image), so chains,
      rings and stars all occur; several sizes can be extracted from one frame.
    * Optional stratification so the set is flat in the radius of gyration or the
      largest centroid separation, with bin edges from equal quantiles of a pilot
      sample or given explicitly, and optional balancing over the contact-graph
      motif.
    * Provenance on every cluster: unique names ``<frame>_<seed>_<molecules>`` and
      ``#ExtractClusters#scan`` properties (size, spread, motif, contacts, bin,
      source molecules) that survive SDF/extxyz export, plus a ``clusters.csv`` per
      step.
