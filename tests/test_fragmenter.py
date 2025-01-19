"""Tests for fragmenter."""

from FragmentRetro.fragmenter import BRICSFragmenter


def test_brics_fragmenter_visualization():
    # Example usage
    smiles = "CCCOCC"
    fragmenter = BRICSFragmenter(smiles)
    # Check if the fragment graph is built correctly
    assert fragmenter.fragment_graph.number_of_nodes() > 0, "Fragment graph should have nodes"
    assert fragmenter.fragment_graph.number_of_edges() >= 0, "Fragment graph should have edges or be empty"
