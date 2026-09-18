# -*- coding: utf-8 -*-
"""
Control parameters for the Extract Clusters step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class ExtractClustersParameters(seamm.Parameters):
    """
    The control parameters for Extract Clusters.

    The keys are the parameters for this plug-in; each value is a dictionary
    describing that parameter (default, kind, units, enumeration, format,
    description, help).

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.

    See Also
    --------
    ExtractClusters, TkExtractClusters, ExtractClustersStep
    """

    parameters = {
        # ------------------------------------------------------------------ #
        # Input: which structures to extract from (standard SEAMM block)
        # ------------------------------------------------------------------ #
        **seamm.standard_parameters.structure_selection_parameters,
        # ------------------------------------------------------------------ #
        # What to extract
        # ------------------------------------------------------------------ #
        "cluster sizes": {
            "default": "3",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Cluster sizes (molecules):",
            "help_text": (
                "The number of molecules per cluster. A comma- or space-separated "
                "list, e.g. '3, 4', extracts a set of each size from the same "
                "frame. Each must be at least 2."
            ),
        },
        "number of clusters": {
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Clusters per size:",
            "help_text": (
                "How many clusters of each size to extract from each selected "
                "configuration (frame). Fewer are produced if the configuration "
                "does not contain enough distinct connected clusters, or if "
                "stratification quotas cannot fill."
            ),
        },
        # ------------------------------------------------------------------ #
        # The contact graph
        # ------------------------------------------------------------------ #
        "contact cutoff": {
            "default": 3.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Contact cutoff:",
            "help_text": (
                "Two molecules are in contact if any of their contact atoms are "
                "within this distance (minimum image). A cluster is a connected "
                "set of molecules under this criterion. 3.5 Å between oxygens is "
                "the usual hydrogen-bond criterion for water."
            ),
        },
        "contact elements": {
            "default": "all",
            "kind": "string",
            "default_units": "",
            "enumeration": ("all",),
            "format_string": "",
            "description": "Contact elements:",
            "help_text": (
                "The elements whose atoms define contact, e.g. 'O' or 'O, N'. "
                "'all' uses every atom. Restricting to the heavy atoms that form "
                "the intermolecular contacts (e.g. O for water) makes the cutoff "
                "meaningful and the graph cheaper to build."
            ),
        },
        # ------------------------------------------------------------------ #
        # Stratification
        # ------------------------------------------------------------------ #
        "spread metric": {
            "default": "radius of gyration",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("radius of gyration", "maximum centroid distance"),
            "format_string": "",
            "description": "Spread metric:",
            "help_text": (
                "The compactness coordinate of a cluster, computed from the "
                "molecular centroids: the radius of gyration, or the largest "
                "centroid-centroid distance. Used for stratification and stored "
                "on each cluster."
            ),
        },
        "stratification": {
            "default": "quantile bins",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("none", "quantile bins", "explicit bin edges"),
            "format_string": "",
            "description": "Stratify by spread:",
            "help_text": (
                "Whether to accept clusters so that the set is flat in the spread "
                "metric. 'quantile bins' places the bin edges at equal quantiles "
                "of a pilot sample of the frame, so the bins adapt to the cluster "
                "size and system; 'explicit bin edges' uses the edges given below; "
                "'none' accepts clusters as sampled (biased towards the compact "
                "clusters that dominate the liquid)."
            ),
        },
        "number of bins": {
            "default": 3,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of bins:",
            "help_text": (
                "The number of equal-quantile bins of the spread metric. The "
                "clusters are accepted so that the bins fill roughly evenly."
            ),
        },
        "bin edges": {
            "default": "",
            "kind": "string",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Bin edges:",
            "help_text": (
                "The bin edges of the spread metric, in Å, as an increasing comma- "
                "or space-separated list, e.g. '1.6, 2.0, 2.5, 3.5'. Clusters "
                "outside the outer edges are rejected. Note that sensible edges "
                "depend strongly on the cluster size and the system; 'quantile "
                "bins' is usually the better choice."
            ),
        },
        "motifs": {
            "default": "any",
            "kind": "string",
            "default_units": "",
            "enumeration": ("any",),
            "format_string": "",
            "description": "Restrict to motifs:",
            "help_text": (
                "Accept only clusters whose contact-graph topology is one of the "
                "given motifs, e.g. 'ring', or 'ring, star' for several. 'any' "
                "accepts all. Trimers have 'chain' and 'ring'; tetramers 'chain', "
                "'star', 'ring', 'paw', 'diamond' and 'K4'; larger clusters are "
                "labelled 'e<n>' by their number of contacts. Rare motifs use up "
                "the attempt budget quickly, so raise 'Attempts per cluster' if "
                "the set comes up short."
            ),
        },
        "balance motifs": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Balance motifs:",
            "help_text": (
                "Also balance the set over the topology of the contact graph "
                "within the cluster (chain / ring / star / ..., or the number of "
                "contacts for clusters of more than four molecules), jointly with "
                "the spread bins. Topology and spread are correlated -- rings only "
                "exist compact -- so some cells cannot fill and the set will "
                "come up short of the requested number; the summary shows which."
            ),
        },
        # ------------------------------------------------------------------ #
        # Sampling
        # ------------------------------------------------------------------ #
        "random seed": {
            "default": "random",
            "kind": "string",
            "default_units": "",
            "enumeration": ("random",),
            "format_string": "",
            "description": "Random seed:",
            "help_text": (
                "The seed for the random-number generator. Use 'random' for a "
                "fresh seed, or an integer for a reproducible selection. The seed "
                "actually used is always printed in the output, so a 'random' run "
                "can be reproduced by entering that value here."
            ),
        },
        "attempts per cluster": {
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Attempts per cluster:",
            "help_text": (
                "The sampling budget: the total number of attempts is this times "
                "the number of clusters requested. Attempts that duplicate an "
                "earlier cluster, fall outside the bin edges, or land in a full "
                "quota cell are rejected and count against the budget."
            ),
        },
        # ------------------------------------------------------------------ #
        # Output
        # ------------------------------------------------------------------ #
        "system name": {
            "default": "from source",
            "kind": "string",
            "default_units": "",
            "enumeration": ("from source",),
            "format_string": "",
            "description": "Name the cluster system:",
            "help_text": (
                "The name of the system that receives the clusters, created if it "
                "does not exist. 'from source' uses '<source system> clusters', "
                "with the system of the first selected structure. The clusters from "
                "all selected frames go into this one system."
            ),
        },
        "name prefix": {
            "default": "from configuration",
            "kind": "string",
            "default_units": "",
            "enumeration": ("from configuration", "none"),
            "format_string": "",
            "description": "Configuration name prefix:",
            "help_text": (
                "Each cluster is named '<prefix><seed>_<m1-m2-...>' from the "
                "source molecules, which is unique within a frame. The prefix "
                "keeps the names unique across frames: 'from configuration' uses "
                "'<source configuration name>_' (the default, and the right choice "
                "when several frames are selected); 'none' uses no prefix; any other "
                "text is used literally (a variable such as '$frame_' works)."
            ),
        },
        "store properties": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Store descriptors as properties:",
            "help_text": (
                "Store the cluster size, spread, motif, number of contacts and "
                "source molecules on each configuration as '#ExtractClusters#scan' "
                "properties. They are carried through SDF/extxyz output."
            ),
        },
        "make current": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Make the cluster system current:",
            "help_text": (
                "Whether to make the cluster system (and its first configuration) "
                "the current one, so that a following Write Structure step writes "
                "the clusters. Choose 'no' when the source must remain current."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ExtractClustersParameters.__init__")

        super().__init__(
            defaults={**ExtractClustersParameters.parameters, **defaults}, data=data
        )
