==========
User Guide
==========

The Extract Clusters step takes the current configuration -- typically a
periodic snapshot of a liquid or solution from Read Structure or a dynamics
step -- and extracts *n*-molecule clusters from it as unwrapped, non-periodic
configurations in a new system. Its main use is generating many-body
(trimer, tetramer, ...) training and diagnostic data for machine-learned
force fields, where a dimer-only training set cannot capture cooperative
interactions.

How clusters are chosen
-----------------------

Molecules are found from the bonds (an ion is a one-atom molecule). Two
molecules are *in contact* if any of their contact atoms -- all atoms, or the
elements you list, e.g. ``O`` for water hydrogen bonds -- are within the
contact cutoff under minimum-image conventions. A cluster is a **connected**
set of *n* molecules under this criterion, grown from a random seed molecule
by random frontier expansion. Every molecule interacts with at least one
other, and because the sampling is random the same frame yields many distinct
clusters with varied topology (chains, rings, stars, ...), unlike selecting
the *n* nearest neighbours, which always gives the most compact cluster.

Stratification
--------------

The compact clusters that dominate a liquid would dominate an unstratified
sample. With *Stratify by spread* set to *quantile bins*, the step draws a
pilot sample of the frame, places bin edges at equal quantiles of the spread
metric (radius of gyration or largest centroid separation of the molecules)
and then accepts clusters so that each bin fills equally. *Explicit bin edges*
lets you give the edges instead; *none* accepts clusters as sampled.

*Balance motifs* also balances the set over the topology of the contact graph
within the cluster. Topology and spread are correlated -- rings only exist
compact -- so some (motif, bin) cells are physically empty and the set then
comes up short of the requested number; the summary table in ``step.out``
shows which cells filled.

Output
------

The clusters go into the system named in *Name the cluster system*
(``<source system> clusters`` by default), each as a non-periodic
configuration centred at the origin with the molecules intact. Names are
``<prefix><seed>_<m1-m2-...>`` from the source molecule indices, so they are
unique within a frame; the prefix (the source configuration name by default)
keeps them unique across frames when a loop extracts from many frames into
one system. The size, spread, motif, number of contacts, bin and source
molecules are stored as ``#ExtractClusters#scan`` properties and written to
``clusters.csv`` in the step directory. A following Write Structure step with
*current system* writes all the clusters, e.g. as an SDF whose records carry
the properties.

Parameters
----------

Cluster sizes (molecules)
    One or more sizes, e.g. ``3`` or ``3, 4``; each at least 2.
Clusters per size
    How many of each size to extract from the frame.
Contact cutoff, Contact elements
    Define the contact graph. 3.5 Å between oxygens is the usual hydrogen-bond
    criterion for water.
Spread metric, Stratify by spread, Number of bins / Bin edges, Balance motifs
    See *Stratification* above. Only the controls that apply to the chosen
    scheme are shown.
Random seed
    ``random`` or an integer for a reproducible selection.
Attempts per cluster
    The sampling budget: total attempts = this × clusters requested.
Name the cluster system, Configuration name prefix, Store descriptors as properties, Make the cluster system current
    See *Output* above. Choose *no* for *Make current* when the source frame
    must stay current, e.g. inside a loop that reads frames into it.
