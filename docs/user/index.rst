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

Which structures
----------------

The step works on the structures chosen in the *Structures to extract from*
block, the standard SEAMM selection: the current configuration (the default),
all or the last or first configurations of the current system, of all systems,
or of the systems whose name is / matches / matches a regular expression; or a
variable holding a list of configurations. Selecting *all* configurations of a
system that holds a trajectory extracts from every frame in one step, so no
loop is needed and all the clusters land in one system, ready for a single
Write Structure. Each frame's clusters are named with that frame's
configuration name as the prefix, and ``clusters.csv`` records the frame.

How clusters are chosen
-----------------------

Molecules are found from the bonds (an ion is a one-atom molecule), so the
structure must carry connectivity: a frame read from an extended XYZ trajectory
needs the Read Structure step's *Perceive bonds* option (the default), and the
step stops with an error if a molecular structure arrives without bonds. Two
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

*Restrict to motifs* accepts only clusters with the given contact-graph topology,
e.g. ``ring`` for cyclic trimers, or several such as ``ring, star``. The dialog
offers the motifs the requested cluster sizes can produce. Rare motifs use up the
attempt budget quickly, so raise *Attempts per cluster* if the set comes up short;
the output reports how many candidates were rejected for their motif.

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
``clusters.csv`` in the step directory; ``summary.json`` there records the seed,
the source, the bin edges and the counts by motif. A following Write Structure step with
*current system* writes all the clusters, e.g. as an SDF whose records carry
the properties.

Parameters
----------

Systems, Configurations
    Which structures to extract from; see *Which structures* above.
Cluster sizes (molecules)
    One or more sizes, e.g. ``3`` or ``3, 4``; each at least 2.
Clusters per size
    How many of each size to extract from each selected frame.
Contact cutoff, Contact elements
    Define the contact graph. 3.5 Å between oxygens is the usual hydrogen-bond
    criterion for water.
Spread metric, Stratify by spread, Number of bins / Bin edges, Balance motifs
    See *Stratification* above. Only the controls that apply to the chosen
    scheme are shown.
Random seed
    ``random`` or an integer. The seed actually used is always printed in the
    output and written to ``summary.json``, so a ``random`` run can be reproduced
    by entering that value.
Attempts per cluster
    The sampling budget: total attempts = this × clusters requested.
Name the cluster system, Configuration name prefix, Store descriptors as properties, Make the cluster system current
    See *Output* above. Choose *no* for *Make current* when the source must
    stay current.
