#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `extract_clusters_step` package."""

import pytest  # noqa: F401
import extract_clusters_step  # noqa: F401
from extract_clusters_step.extract_clusters import ExtractClusters


def test_construction():
    """Just create an object and test its type."""
    result = extract_clusters_step.ExtractClusters()
    assert (
        str(type(result))
        == "<class 'extract_clusters_step.extract_clusters.ExtractClusters'>"
    )


def test_description_text():
    node = extract_clusters_step.ExtractClusters()
    node._id = (1,)  # normally assigned by the flowchart; needed for the header
    text = " ".join(node.description_text().split())  # undo the line wrapping
    assert "Extract 50 clusters of 3 molecule" in text
    assert "equal-quantile bins" in text
    P = node.parameters.values_to_dict()
    P["cluster sizes"] = "3, 4"
    P["stratification"] = "none"
    P["balance motifs"] = "yes"
    P["random seed"] = "42"
    P["system name"] = "tets"
    P["name prefix"] = "none"
    text = " ".join(node.description_text(P).split())  # undo the line wrapping
    assert "3, 4 molecules" in text
    assert "without stratification" in text
    assert "balanced over the topology" in text
    assert "random seed is 42" in text
    assert "'tets'" in text


def test_parse_sizes():
    assert ExtractClusters._parse_sizes("3") == [3]
    assert ExtractClusters._parse_sizes("3, 4 5,3") == [3, 4, 5]
    assert ExtractClusters._parse_sizes([3, "4"]) == [3, 4]
    with pytest.raises(ValueError):
        ExtractClusters._parse_sizes("1")
    with pytest.raises(ValueError):
        ExtractClusters._parse_sizes("three")
    with pytest.raises(ValueError):
        ExtractClusters._parse_sizes("")


def test_parse_elements():
    assert ExtractClusters._parse_elements("all") is None
    assert ExtractClusters._parse_elements("") is None
    assert ExtractClusters._parse_elements("O") == ["O"]
    assert ExtractClusters._parse_elements("o, n Li") == ["O", "N", "Li"]
    with pytest.raises(ValueError):
        ExtractClusters._parse_elements("Xx")


def test_parse_edges():
    assert ExtractClusters._parse_edges("1.6, 2.0 2.5") == [1.6, 2.0, 2.5]
    with pytest.raises(ValueError):
        ExtractClusters._parse_edges("1.6")
    with pytest.raises(ValueError):
        ExtractClusters._parse_edges("2.0, 1.6")
    with pytest.raises(ValueError):
        ExtractClusters._parse_edges("a, b")


def test_make_rng():
    import numpy as np

    r1 = ExtractClusters._make_rng("42").integers(1000)
    r2 = ExtractClusters._make_rng(42).integers(1000)
    assert r1 == r2
    assert isinstance(ExtractClusters._make_rng("random"), np.random.Generator)
