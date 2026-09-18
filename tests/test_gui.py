# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice
that drives the layout. Skipped when no display is available."""

import pytest


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Extract Clusters")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Extract Clusters")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    for strat in ("none", "explicit bin edges", "quantile bins"):
        tk_node["stratification"].set(strat)
        tk_node.reset_dialog()
    # sizes drive the motif list; a single motif hides 'balance motifs'
    tk_node["cluster sizes"].set("4")
    tk_node.reset_dialog()
    assert "K4" in tk_node["motifs"].combobox.cget("values")
    tk_node["motifs"].set("ring")
    tk_node.reset_dialog()
    assert not tk_node["balance motifs"].winfo_ismapped()
    tk_node["motifs"].set("any")
    tk_node.reset_dialog()
    for systems in ("name matches", "all", "current"):
        tk_node["source systems"].set(systems)
        tk_node.reset_dialog()
    # the name field appears only for a by-name choice
    tk_node["source configurations"].set("name is")
    tk_node.reset_dialog()
    assert tk_node["source configuration name"].grid_info() != {}
    tk_node["source configurations"].set("all")
    tk_node.reset_dialog()
    assert tk_node["source configuration name"].grid_info() == {}
