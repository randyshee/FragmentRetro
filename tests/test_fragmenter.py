"""Tests for fragmenter."""

from itertools import combinations as itertools_combinations

import pytest

from FragmentRetro.fragmenter import BRICSFragmenter

TEST_CASES_FOR_GET_LENGTH_N_COMBINATIONS = [
    {
        "case_number": 1,
        "smiles": "CNCc1cc(-c2ccccc2F)n(S(=O)(=O)c2cccnc2)c1",
    },
    {
        "case_number": 2,
        "smiles": "CCCCCN(CCCC)CCCN(CCOC)c1ccccc1",
    },
    {
        "case_number": 3,
        "smiles": "COc1ccc(-n2nccn2)c(C(=O)N2CCC[C@@]2(C)c2nc3c(C)c(Cl)ccc3[nH]2)c1",
    },
]


def test_brics_fragmenter_visualization():
    smiles = "CCCOCC"
    fragmenter = BRICSFragmenter(smiles)
    # Check if the fragment graph is built correctly
    assert fragmenter.fragment_graph.number_of_nodes() > 0, "Fragment graph should have nodes"
    assert fragmenter.fragment_graph.number_of_edges() >= 0, "Fragment graph should have edges or be empty"


@pytest.mark.parametrize(
    "case_number, smiles",
    [
        (
            tc["case_number"],
            tc["smiles"],
        )
        for tc in TEST_CASES_FOR_GET_LENGTH_N_COMBINATIONS
    ],
)
def test_get_length_n_combinations(case_number, smiles):
    fragmenter = BRICSFragmenter(smiles)
    num_fragments = fragmenter.num_fragments
    for length_n in range(1, num_fragments + 1):
        good_combinations = fragmenter.get_length_n_combinations(length_n)
        all_combinations = list(itertools_combinations(range(num_fragments), length_n))
        bad_combinations = [comb for comb in all_combinations if comb not in good_combinations]
        for good_comb in good_combinations:
            assert fragmenter.check_connected_subgraph(
                good_comb
            ), f"Case {case_number} Length {length_n}: Good combination {good_comb} should be connected"
            assert (
                len(set(good_comb)) == length_n
            ), f"Case {case_number} Length {length_n}: Good combination {good_comb} should have repeated elements"
        for bad_comb in bad_combinations:
            assert not fragmenter.check_connected_subgraph(
                bad_comb
            ), f"Case {case_number} Length {length_n}: Bad combination {bad_comb} should not be connected"


if __name__ == "__main__":
    pytest.main()
