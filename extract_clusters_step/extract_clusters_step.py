# -*- coding: utf-8 -*-

import extract_clusters_step


class ExtractClustersStep(object):
    """Helper class needed for the stevedore integration.

    This must provide a description() method that returns a dict containing a
    description of this node, and create_node() and create_tk_node() methods
    for creating the graphical and non-graphical nodes.

    Attributes
    ----------
    my_description : {description, group, name}
        A human-readable description of this step. It can be
        several lines long, and needs to be clear to non-expert users.
        It contains the following keys: description, group, name.
    my_description["description"] : tuple
        A description of the Extract Clusters step. It must be
        clear to non-experts.
    my_description["group"] : str
        Which group in the menus to put this step. If the group does
        not exist it will be created. Common groups are "Building",
        "Calculations", "Control" and "Data".
    my_description["name"] : str
        The name of this step, to be displayed in the menus.
    """

    my_description = {
        "description": (
            "Extract molecular clusters (trimers, tetramers, ... larger n-mers) "
            "from a condensed-phase, typically periodic, configuration as "
            "unwrapped, non-periodic structures, e.g. for many-body training data "
            "and diagnostics for machine-learned force fields."
        ),
        "group": "Building",
        "name": "Extract Clusters",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does."""
        return ExtractClustersStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Create and return the new node object.

        Parameters
        ----------
        flowchart: seamm.Node
            A non-graphical SEAMM node

        **kwargs : keyworded arguments
            Various keyworded arguments such as title, namespace or
            extension representing the title displayed in the flowchart,
            the namespace for the plugins of a subflowchart and
            the extension, respectively.

        Returns
        -------
        ExtractClusters

        See Also
        --------
        ExtractClusters

        """

        return extract_clusters_step.ExtractClusters(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Create and return the graphical Tk node object.

        Parameters
        ----------
        canvas : tk.Canvas
            The Tk Canvas widget

        **kwargs : keyworded arguments
            Various keyworded arguments such as tk_flowchart, node, x, y, w, h
            representing a graphical flowchart object, a non-graphical node for
            a step, and dimensions of the graphical node.

        Returns
        -------
        TkExtractClusters

        See Also
        --------
        TkExtractClusters
        """

        return extract_clusters_step.TkExtractClusters(canvas=canvas, **kwargs)
