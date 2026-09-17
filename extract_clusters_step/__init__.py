# -*- coding: utf-8 -*-

"""
extract_clusters_step
A SEAMM plug-in for extracting molecular clusters from a periodic cell of molecules
"""

# Bring up the classes so that they appear to be directly in
# the extract_clusters_step package.

from extract_clusters_step.extract_clusters import ExtractClusters  # noqa: F401, E501
from extract_clusters_step.extract_clusters_parameters import ExtractClustersParameters  # noqa: F401, E501
from extract_clusters_step.extract_clusters_step import ExtractClustersStep  # noqa: F401, E501
from extract_clusters_step.tk_extract_clusters import TkExtractClusters  # noqa: F401, E501

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
