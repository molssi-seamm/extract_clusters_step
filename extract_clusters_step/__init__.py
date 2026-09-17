# -*- coding: utf-8 -*-

"""
extract_clusters_step
A SEAMM plug-in for extracting molecular clusters from a periodic cell of molecules
"""

# Bring up the classes so that they appear to be directly in
# the extract_clusters_step package.

from extract_clusters_step.extract_clusters import ExtractClusters  # noqa: F401, E501
from extract_clusters_step.extract_clusters_parameters import (  # noqa: F401
    ExtractClustersParameters,
)
from extract_clusters_step.extract_clusters_step import (  # noqa: F401
    ExtractClustersStep,
)
from extract_clusters_step.tk_extract_clusters import (  # noqa: F401
    TkExtractClusters,
)
from extract_clusters_step.metadata import metadata  # noqa: F401
from extract_clusters_step.cluster_sampling import (  # noqa: F401
    PROPERTY_TAG,
    classify_motif,
    cluster_summary,
    extract_nmers,
)

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
