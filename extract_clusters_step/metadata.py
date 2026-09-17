# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Extract Clusters"""

metadata = {}

# Extract Clusters generates structures rather than computing properties, so it
# exposes no formal results. The per-configuration descriptors (size, spread,
# motif, source molecules) are stored directly as properties tagged
# '#ExtractClusters#scan' -- see data/properties.csv. An empty dict is still
# required by the framework's results handling.
metadata["results"] = {}
