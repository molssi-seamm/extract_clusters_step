#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `extract_clusters_step` package."""

import pytest  # noqa: F401
import extract_clusters_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = extract_clusters_step.ExtractClusters()
    assert str(type(result)) == "<class 'extract_clusters_step.extract_clusters.ExtractClusters'>"
