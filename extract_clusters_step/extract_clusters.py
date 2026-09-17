# -*- coding: utf-8 -*-

"""Non-graphical part of the Extract Clusters step in a SEAMM flowchart"""

import csv
import importlib.resources
import logging
from pathlib import Path
import pprint  # noqa: F401
import re

import numpy as np

import extract_clusters_step
from extract_clusters_step.cluster_sampling import (
    PROPERTY_TAG,
    cluster_summary,
    extract_nmers,
)
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Extract Clusters")

# Add this module's properties to the standard properties
path = importlib.resources.files("extract_clusters_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

_spread_metric = {
    "radius of gyration": "rg",
    "maximum centroid distance": "dmax",
}


class ExtractClusters(seamm.Node):
    """
    The non-graphical part of a Extract Clusters step in a flowchart.

    Extracts n-molecule clusters (trimers, tetramers, ... larger clusters) from
    the current, typically periodic, condensed-phase configuration as unwrapped,
    non-periodic configurations in a new system, tagged with provenance. The
    clusters are connected subgraphs of the molecular contact graph, optionally
    stratified to be flat in a spread coordinate and balanced over topology.
    See :mod:`extract_clusters_step.cluster_sampling` for the algorithm.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    parameters : ExtractClustersParameters
        The control parameters for Extract Clusters.

    See Also
    --------
    TkExtractClusters, ExtractClustersParameters
    """

    def __init__(
        self, flowchart=None, title="Extract Clusters", extension=None, logger=logger
    ):
        """A step for Extract Clusters in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Extract Clusters {self}")

        super().__init__(
            flowchart=flowchart,
            title="Extract Clusters",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self.parameters = extract_clusters_step.ExtractClustersParameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return extract_clusters_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return extract_clusters_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        sizes = str(P["cluster sizes"]).strip()
        plural = "s" if ("," in sizes or " " in sizes) else ""
        text = (
            f"Extract {P['number of clusters']} clusters of {sizes} molecule{plural}"
            " each from the current configuration. Molecules are in contact if "
            f"their contact atoms ({P['contact elements']}) are within "
            f"{P['contact cutoff']}, and a cluster is a connected set of "
            "molecules under that criterion, sampled by random growth from a "
            "random seed molecule."
        )

        strat = P["stratification"]
        metric = P["spread metric"]
        if strat == "quantile bins":
            text += (
                f" The clusters will be stratified to be flat in the {metric}, "
                f"using {P['number of bins']} equal-quantile bins from a pilot "
                "sample of the frame."
            )
        elif strat == "explicit bin edges":
            text += (
                f" The clusters will be stratified to be flat in the {metric}, "
                f"using the bin edges {P['bin edges']} Å."
            )
        else:
            text += " The clusters are accepted as sampled, without stratification."
        if self._truthy(P["balance motifs"]):
            text += (
                " The set will also be balanced over the topology of the contact "
                "graph (chain, ring, star, ...)."
            )

        seed = str(P["random seed"]).strip()
        if seed.lower() in ("", "random"):
            text += " The random seed will be chosen at random."
        else:
            text += f" The random seed is {seed}."

        if P["system name"] == "from source":
            text += (
                " The clusters will be placed as unwrapped, non-periodic "
                "configurations in the system '<source system> clusters',"
            )
        else:
            text += (
                " The clusters will be placed as unwrapped, non-periodic "
                f"configurations in the system '{P['system name']}',"
            )
        if P["name prefix"] == "from configuration":
            text += " named '<source configuration>_<seed>_<molecules>'."
        elif P["name prefix"] == "none":
            text += " named '<seed>_<molecules>'."
        else:
            text += f" named '{P['name prefix']}<seed>_<molecules>'."
        if self._truthy(P["store properties"]):
            text += (
                " The size, spread, motif and source molecules of each cluster "
                "will be stored as properties on its configuration."
            )
        if self._truthy(P["make current"]):
            text += " The cluster system will be made the current system."

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Extract Clusters step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Print what we are doing
        printer.important(__(self.description_text(P), indent=self.indent))
        printer.important("")

        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # The source: the current system and configuration
        system_db = self.get_variable("_system_db")
        system, configuration = self.get_system_configuration(None)
        if configuration is None or configuration.n_atoms == 0:
            raise RuntimeError("Extract Clusters: there is no current structure.")

        # Parse the control parameters (run-time backstop for the GUI checks)
        sizes = self._parse_sizes(P["cluster sizes"])
        n_samples = int(P["number of clusters"])
        if n_samples < 1:
            raise ValueError("The number of clusters must be at least 1.")
        cutoff = P["contact cutoff"]
        if isinstance(cutoff, Q_):
            cutoff = cutoff.m_as("Å")
        cutoff = float(cutoff)
        if cutoff <= 0:
            raise ValueError("The contact cutoff must be positive.")
        contact_elements = self._parse_elements(P["contact elements"])
        spread_metric = _spread_metric[P["spread metric"]]
        strat = P["stratification"]
        if strat == "quantile bins":
            spread_bins = int(P["number of bins"])
            if spread_bins < 1:
                raise ValueError("The number of bins must be at least 1.")
        elif strat == "explicit bin edges":
            spread_bins = self._parse_edges(P["bin edges"])
        elif strat == "none":
            spread_bins = None
        else:
            raise ValueError(f"Unknown stratification '{strat}'")
        balance_motifs = self._truthy(P["balance motifs"])
        rng = self._make_rng(P["random seed"])
        attempts = int(P["attempts per cluster"])
        if attempts < 1:
            raise ValueError("The attempts per cluster must be at least 1.")
        store_properties = self._truthy(P["store properties"])

        # Destination system and naming
        system_name = P["system name"]
        if system_name == "from source":
            system_name = f"{system.name} clusters"
        systems = system_db.get_systems(system_name)
        if systems:
            new_system = systems[0]
        else:
            new_system = system_db.create_system(name=system_name, make_current=False)

        prefix = str(P["name prefix"])
        if prefix == "from configuration":
            prefix = f"{configuration.name}_"
        elif prefix == "none":
            prefix = ""

        # Extract each size, sharing the 'seen' set so nothing is duplicated
        results = []
        seen = set()
        for n in sizes:
            configurations, records, info = extract_nmers(
                configuration,
                n,
                n_samples,
                cutoff=cutoff,
                contact_elements=contact_elements,
                spread_metric=spread_metric,
                spread_bins=spread_bins,
                balance_motifs=balance_motifs,
                system=new_system,
                name_prefix=prefix,
                max_attempts=attempts * n_samples,
                rng=rng,
                store_properties=store_properties,
                seen=seen,
            )
            results.append((n, configurations, records, info))

        # Record the descriptors for later analysis
        self._write_csv(directory / "clusters.csv", results)

        # Make the new system & its first configuration current
        if self._truthy(P["make current"]) and new_system.n_configurations > 0:
            system_db.system = new_system
            new_system.configuration = new_system.configurations[0].id

        self.analyze(P=P, results=results, system=new_system, source=configuration)

        return next_node

    def analyze(
        self, indent="", P=None, results=None, system=None, source=None, **kwargs
    ):
        """Report the extracted clusters to step.out.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        P : dict
            The control parameters.
        results : [(n, configurations, records, info)]
            The output of :func:`extract_nmers` for each cluster size.
        system : molsystem _System
            The destination system.
        source : molsystem _Configuration
            The source configuration.
        """
        if results is None:
            return

        metric = P["spread metric"] if P is not None else "spread"
        total = 0
        for n, configurations, records, info in results:
            total += len(configurations)
            text = (
                f"{n}-mers: extracted {len(configurations)} clusters in "
                f"{info['attempts']} attempts from {info['n_molecules']} molecules "
                f"with {info['n_contacts']} contacts."
            )
            if records:
                rg = [r["rg"] for r in records]
                dmax = [r["dmax"] for r in records]
                text += (
                    f" Radius of gyration {min(rg):.2f}-{max(rg):.2f} Å; largest "
                    f"centroid separation {min(dmax):.2f}-{max(dmax):.2f} Å."
                )
            printer.important(__(text, indent=4 * " "))
            if records:
                printer.important("")
                printer.important(
                    __(
                        f"Clusters by motif and {metric} bin:",
                        indent=8 * " ",
                        wrap=False,
                    )
                )
                table = cluster_summary(records, edges=info["edges"])
                for line in table.splitlines():
                    printer.important(
                        __(line, indent=12 * " ", wrap=False, dedent=False)
                    )
            if info["warning"] is not None:
                printer.important("")
                printer.important(
                    __("Warning: " + info["warning"], indent=8 * " ", wrap=True)
                )
            printer.important("")

        if system is not None:
            text = (
                f"In total {total} clusters were added to the system "
                f"'{system.name}', which now has {system.n_configurations} "
                "configurations."
            )
            if source is not None:
                text += (
                    f" The source was configuration '{source.name}' of system "
                    f"'{source.system.name}'."
                )
            printer.important(__(text, indent=4 * " "))
        printer.important(
            __(
                "The descriptors of each cluster are in 'clusters.csv' in this "
                "step's directory.",
                indent=4 * " ",
            )
        )
        printer.important("")

    # ----------------------------------------------------------------- #
    # Implementation helpers
    # ----------------------------------------------------------------- #

    @staticmethod
    def _split_list(text):
        """Split a comma- and/or whitespace-separated string into tokens."""
        if isinstance(text, (list, tuple)):
            return [str(t).strip() for t in text if str(t).strip() != ""]
        return [t for t in re.split(r"[,\s]+", str(text).strip()) if t != ""]

    @classmethod
    def _parse_sizes(cls, text):
        """The list of cluster sizes from the 'cluster sizes' parameter."""
        tokens = cls._split_list(text)
        if not tokens:
            raise ValueError("At least one cluster size must be given.")
        sizes = []
        for t in tokens:
            try:
                n = int(t)
            except ValueError:
                raise ValueError(f"The cluster size '{t}' is not an integer.")
            if n < 2:
                raise ValueError(f"The cluster size must be at least 2, not {n}.")
            if n not in sizes:
                sizes.append(n)
        return sizes

    @classmethod
    def _parse_elements(cls, text):
        """The contact elements, or None for all atoms."""
        if isinstance(text, str) and text.strip().lower() in ("", "all"):
            return None
        tokens = cls._split_list(text)
        if not tokens:
            return None
        elements = []
        for t in tokens:
            symbol = t.capitalize()
            if symbol not in molsystem.elements.symbol_to_atno:
                raise ValueError(f"'{t}' is not an element symbol.")
            if symbol not in elements:
                elements.append(symbol)
        return elements

    @classmethod
    def _parse_edges(cls, text):
        """The explicit bin edges (Å) as an increasing list of floats."""
        tokens = cls._split_list(text)
        if len(tokens) < 2:
            raise ValueError(
                "Explicit stratification needs at least two bin edges, e.g. "
                "'1.6, 2.0, 2.5'."
            )
        edges = []
        for t in tokens:
            try:
                edges.append(float(t))
            except ValueError:
                raise ValueError(f"The bin edge '{t}' is not a number.")
        if any(b <= a for a, b in zip(edges, edges[1:])):
            raise ValueError("The bin edges must be strictly increasing.")
        return edges

    @staticmethod
    def _truthy(value):
        """A boolean parameter is a bool at run time but 'yes'/'no' in the GUI."""
        return value is True or (isinstance(value, str) and value.lower() == "yes")

    @staticmethod
    def _make_rng(seed):
        """A numpy random generator from the 'random seed' parameter."""
        if isinstance(seed, str):
            if seed.strip() == "" or seed.strip().lower() == "random":
                return np.random.default_rng()
            seed = int(seed)
        return np.random.default_rng(int(seed))

    @staticmethod
    def _write_csv(path, results):
        """Write one row per cluster with its descriptors."""
        fields = [
            "name",
            "size",
            "motif",
            "bin",
            "rg",
            "dmax",
            "n_edges",
            "degrees",
            "seed",
            "molecules",
        ]
        with open(path, "w", newline="") as fd:
            writer = csv.writer(fd)
            writer.writerow(fields)
            for n, configurations, records, info in results:
                for r in records:
                    writer.writerow(
                        [
                            r["name"],
                            r["n"],
                            r["motif"],
                            r["bin"],
                            f"{r['rg']:.4f}",
                            f"{r['dmax']:.4f}",
                            r["n_edges"],
                            "-".join(str(d) for d in r["degrees"]),
                            r["seed"],
                            "-".join(str(m) for m in r["molecules"]),
                        ]
                    )

    @property
    def property_tag(self):
        """The suffix of the per-cluster properties, '#ExtractClusters#scan'."""
        return PROPERTY_TAG
