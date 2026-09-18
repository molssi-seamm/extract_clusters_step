=========================================================
2026-09-18: Reproducibility, motif selection, many frames
=========================================================

Three enhancements requested after the first real use of the step on the
3333-water NVT trajectory (Job 4800), planned together and delivered in two
releases.

1. Print the random seed
------------------------

A run with the seed set to ``random`` could not be reproduced. The step now
draws the seed from the operating system, builds the generator from it, prints
it in the output and records it in ``summary.json`` together with the source
system and configuration, the bin edges and the counts by motif. Entering the
printed value as the seed reproduces the selection exactly (same source, same
parameters). Delivered in 2026.9.18.

2. Restrict the topology
------------------------

A new *Restrict to motifs* parameter accepts only clusters whose contact-graph
motif is in the given list, e.g. ``ring`` or ``ring, star``. The names come from
``motif_names(n)`` in ``cluster_sampling``, so the dialog can offer exactly the
motifs the requested sizes can produce; the run-time check validates typed or
scripted values the same way. Rejection happens before the stratification
quotas, uses the attempt budget, and is counted and reported, since rare motifs
(rings are a few percent of water trimers) otherwise look like a sampling
failure. A single motif hides *Balance motifs*, which would have nothing to
balance over. Delivered in 2026.9.18.

3. Operate on many configurations in one step
---------------------------------------------

The step is almost always used in a loop over trajectory frames, which is
painful for the user and fragments the output into one directory and one file
per iteration. The same selection problem has been solved three times already
with three vocabularies: the Loop step (system-name filter plus "default
configuration"), Write Structure ("current configuration / current system /
all systems" plus a configuration filter) and the Dimer Builder
(``_resolve_pool``: current, a system name, or a ``$variable`` plus
all / last / first / name is / matches / regexp).

Part A -- a shared selection facility in ``seamm.Node``, next to
``get_system_configuration``, since every step inherits it and it owns the
variables and the current system:

* a standard parameter block (which structures: current configuration, current
  system, a system, all systems, a variable; a system-name filter as in Loop; a
  configuration filter as in Loop and the Dimer Builder) that a step includes in
  its own parameters like the structure-handling block;
* ``Node.select_configurations(P)`` returning the ordered list, resolving
  ``$variables`` holding lists of configurations, with clear errors for an
  empty selection;
* a ``TkNode`` helper that lays the block out reactively so every dialog looks
  the same;
* then convert Loop (the reference semantics), Write Structure and the Dimer
  Builder one at a time, keeping flowchart parameter names backward compatible.

Part A was delivered in seamm 2026.9.18.1 (``structure_selection_parameters``,
``select_configurations``, ``structure_selection_description``,
``TkNode.create_structure_selection_widgets`` / ``layout_structure_selection``,
developer-guide page "Selecting structures in a plug-in"). While reading the
Loop step's selection code two bugs were noted for its conversion: it compares
configuration names with ``is`` rather than ``==``, and it tests for
``matches`` / ``regexp`` while its enumeration offers ``name matches`` /
``name regexp``, so those choices fall through.

Part B -- Extract Clusters uses the block (delivered in 2026.9.18.1). Per selected configuration the
extraction runs with the configuration name as the prefix, sharing one
destination system and one random stream so the printed seed reproduces the
whole set; the dedupe set stays per frame (the same molecule indices in
different frames are different geometries); ``clusters.csv`` gains a ``frame``
column; the report gives one motif-by-bin table per size aggregated over frames
plus a per-frame line.

Verified on six frames of the Job 4800 trajectory read into one system: 60
trimers and 60 tetramers in one step, 10 of each per frame, one SDF.

Sequencing: Part A is a ``seamm`` release with tests against an in-memory
database and a write-up of the selection semantics in the seamm developer
docs; Part B follows as an Extract Clusters release once Part A is on PyPI;
the conversions of Loop and Write Structure come after, each its own release.

Open decisions: whether Part A lives in ``seamm.Node`` (preferred: it needs the
variables and the current system) or a separate library; whether the motif
restriction takes a list (yes -- implemented in item 2).
