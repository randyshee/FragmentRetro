"""Tests for fragmenter."""

from itertools import combinations as itertools_combinations

import pytest

from FragmentRetro.fragmenter import BRICSFragmenter, rBRICSFragmenter
from FragmentRetro.substructure_matcher import SubstructureMatcher
from FragmentRetro.utils.helpers import count_heavy_atoms, replace_dummy_atoms_regex

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


@pytest.mark.parametrize("fragmenter_class", [BRICSFragmenter, rBRICSFragmenter])
def test_fragmenter_visualization(fragmenter_class):
    smiles = "CCCOCC"
    fragmenter = fragmenter_class(smiles)
    # Check if the fragment graph is built correctly
    assert fragmenter.fragment_graph.number_of_nodes() > 0, "Fragment graph should have nodes"
    assert fragmenter.fragment_graph.number_of_edges() >= 0, "Fragment graph should have edges or be empty"


@pytest.mark.parametrize(
    "fragmenter_class, case_number, smiles",
    [(BRICSFragmenter, tc["case_number"], tc["smiles"]) for tc in TEST_CASES_FOR_GET_LENGTH_N_COMBINATIONS]
    + [(rBRICSFragmenter, tc["case_number"], tc["smiles"]) for tc in TEST_CASES_FOR_GET_LENGTH_N_COMBINATIONS],
)
def test_get_length_n_combinations(fragmenter_class, case_number, smiles):
    fragmenter = fragmenter_class(smiles)
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


@pytest.mark.parametrize(
    "fragmenter_class, case_number, smiles",
    [
        (fragmenter_class, tc["case_number"], tc["smiles"])
        for fragmenter_class in [BRICSFragmenter, rBRICSFragmenter]
        for tc in TEST_CASES_FOR_GET_LENGTH_N_COMBINATIONS
    ],
)
def test_get_combination_smiles(fragmenter_class, case_number, smiles):
    fragmenter = fragmenter_class(smiles)
    num_fragments = fragmenter.num_fragments
    for length_n in range(2, num_fragments + 1):
        large_good_combinations = fragmenter.get_length_n_combinations(length_n)
        small_good_combinations = fragmenter.get_length_n_combinations(length_n - 1)
        for large_good_comb in large_good_combinations:
            large_fragment_smiles = fragmenter.get_combination_smiles(large_good_comb)
            large_no_dummy = replace_dummy_atoms_regex(large_fragment_smiles)
            for small_good_comb in small_good_combinations:
                if set(small_good_comb).issubset(set(large_good_comb)):
                    small_fragment_smiles = fragmenter.get_combination_smiles(small_good_comb)
                    small_no_dummy = replace_dummy_atoms_regex(small_fragment_smiles)
                    assert (
                        count_heavy_atoms(small_no_dummy) < count_heavy_atoms(large_no_dummy)
                    ), f"Case {case_number} Length {length_n}: small fragment {small_fragment_smiles} should have fewer heavy atoms than large fragment {large_fragment_smiles}"
                    assert SubstructureMatcher.is_strict_substructure(
                        small_fragment_smiles, large_fragment_smiles
                    ), f"Case {case_number} Length {length_n}: small fragment {small_fragment_smiles} should be a strict substructure of large fragment {large_fragment_smiles}"


if __name__ == "__main__":
    pytest.main()
