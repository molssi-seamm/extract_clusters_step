# -*- coding: utf-8 -*-

"""The graphical part of a Extract Clusters step"""

import pprint  # noqa: F401
import re
import tkinter as tk
import tkinter.ttk as ttk

import extract_clusters_step  # noqa: F401
from extract_clusters_step.cluster_sampling import motif_names
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkExtractClusters(seamm.TkNode):
    """
    The graphical part of a Extract Clusters step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Extract Clusters_parameters.py

    See Also
    --------
    ExtractClusters, TkExtractClusters,
    ExtractClustersParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Extract Clusters_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkExtractClusters.reset_dialog
        """

        frame = super().create_dialog(title="Extract Clusters")
        # Shortcut for parameters
        P = self.node.parameters

        # Group the widgets into labeled frames
        self["clusters frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Clusters",
            labelanchor="n",
            padding=10,
        )
        self["stratification frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Stratification",
            labelanchor="n",
            padding=10,
        )
        self["output frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Output",
            labelanchor="n",
            padding=10,
        )
        parents = {
            "cluster sizes": "clusters frame",
            "number of clusters": "clusters frame",
            "contact cutoff": "clusters frame",
            "contact elements": "clusters frame",
            "random seed": "clusters frame",
            "attempts per cluster": "clusters frame",
            "spread metric": "stratification frame",
            "stratification": "stratification frame",
            "number of bins": "stratification frame",
            "bin edges": "stratification frame",
            "motifs": "stratification frame",
            "balance motifs": "stratification frame",
            "system name": "output frame",
            "name prefix": "output frame",
            "store properties": "output frame",
            "make current": "output frame",
        }

        # Then create the widgets
        for key in P:
            self[key] = P[key].widget(self[parents[key]])

        # Comboboxes whose value changes the layout re-lay out the dialog.
        for key in ("stratification", "motifs"):
            self[key].combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].combobox.bind("<Return>", self.reset_dialog)
            self[key].combobox.bind("<FocusOut>", self.reset_dialog)
        # The cluster sizes decide which motifs exist.
        self["cluster sizes"].entry.bind("<Return>", self.reset_dialog)
        self["cluster sizes"].entry.bind("<FocusOut>", self.reset_dialog)

        # and lay them out
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        Only the stratification controls that apply to the chosen scheme are
        shown: the number of bins for 'quantile bins', the edges for 'explicit
        bin edges', and neither when stratification is off, so an unusable
        combination cannot be built in the editor.

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkExtractClusters.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()
        for name in ("clusters frame", "stratification frame", "output frame"):
            for slave in self[name].grid_slaves():
                slave.grid_forget()

        # keep track of the row in a variable, so that the layout is flexible
        # if e.g. rows are skipped to control such as "method" here
        row = 0
        self["clusters frame"].grid(row=row, column=0, sticky=tk.EW, pady=5)
        row += 1
        self["stratification frame"].grid(row=row, column=0, sticky=tk.EW, pady=5)
        row += 1
        self["output frame"].grid(row=row, column=0, sticky=tk.EW, pady=5)
        row += 1
        frame.columnconfigure(0, weight=1)

        # The clusters
        widgets = []
        for i, key in enumerate(
            (
                "cluster sizes",
                "number of clusters",
                "contact cutoff",
                "contact elements",
                "random seed",
                "attempts per cluster",
            )
        ):
            self[key].grid(row=i, column=0, sticky=tk.EW)
            widgets.append(self[key])
        sw.align_labels(widgets, sticky=tk.E)

        # Stratification, showing only the controls the scheme needs
        strat = self["stratification"].get()
        keys = ["spread metric", "stratification"]
        if strat == "quantile bins":
            keys.append("number of bins")
        elif strat == "explicit bin edges":
            keys.append("bin edges")

        # Offer only the motifs the requested cluster sizes can produce (a
        # typed $variable or unparsable sizes fall back to the n = 3, 4 names).
        sizes = []
        for token in re.split(r"[,\s]+", self["cluster sizes"].get().strip()):
            if token.isdigit() and int(token) >= 2:
                sizes.append(int(token))
        if not sizes:
            sizes = [3, 4]
        names = ["any"]
        for n in sizes:
            names += [m for m in motif_names(n) if m not in names]
        self["motifs"].combobox.config(values=names)
        motifs = [t for t in re.split(r"[,\s]+", self["motifs"].get().strip()) if t]
        keys.append("motifs")
        # A single motif leaves nothing to balance over, so hide the option.
        if not (len(motifs) == 1 and motifs[0].lower() != "any"):
            keys.append("balance motifs")
        widgets = []
        for i, key in enumerate(keys):
            self[key].grid(row=i, column=0, sticky=tk.EW)
            widgets.append(self[key])
        sw.align_labels(widgets, sticky=tk.E)

        # Output
        widgets = []
        for i, key in enumerate(
            ("system name", "name prefix", "store properties", "make current")
        ):
            self[key].grid(row=i, column=0, sticky=tk.EW)
            widgets.append(self[key])
        sw.align_labels(widgets, sticky=tk.E)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkExtractClusters.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def edit(self):
        """Present a dialog for editing the Extract Clusters input

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkExtractClusters.right_click
        """

        if self.dialog is None:
            self.create_dialog()

        self.dialog.activate(geometry="centerscreenfirst")

    def handle_dialog(self, result):
        """Handle the closing of the edit dialog

        What to do depends on the button used to close the dialog. If
        the user closes it by clicking the "x" of the dialog window,
        None is returned, which we take as equivalent to cancel.

        Parameters
        ----------
        result : None or str
            The value of this variable depends on what the button
            the user clicked.

        Returns
        -------
        None
        """

        if result is None or result == "Cancel":
            self.dialog.deactivate(result)
            return

        if result == "Help":
            # display help!!!
            return

        if result != "OK":
            self.dialog.deactivate(result)
            raise RuntimeError(f"Don't recognize dialog result '{result}'")

        self.dialog.deactivate(result)
        # Shortcut for parameters
        P = self.node.parameters

        # Get the values for all the widgets. This may be overkill, but
        # it is easy! You can sort out what it all means later, or
        # be a bit more selective.
        for key in P:
            P[key].set_from_widget()

    def handle_help(self):
        """Shows the help to the user when click on help button.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Help not implemented yet for Extract Clusters!")
