"""Tests for compound filter."""

from pathlib import Path

import pytest

from FragmentRetro.substructure_matcher import SubstructureMatcher
from FragmentRetro.utils.filter_compound import CompoundFilter, precompute_properties

DATA_PATH = Path(__file__).parent.parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
PRECOMPUTE_PATH = DATA_PATH / "precompute"
MOL_PROPERTIES_PATH = PRECOMPUTE_PATH / "n1_stock_properties_subset.json"


with open(PAROUTES_PATH / "n1-stock.txt", "r") as f:
    n1_stock = [line.strip() for line in f.readlines()]
n1_stock_subset = n1_stock[:500]
PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
precompute_properties(n1_stock_subset, MOL_PROPERTIES_PATH, fpSize=2048)

fragment_smiles_list = ["[5*]N1CCC[C@@]1([13*])C", "[4*]CCN[5*]", "[4*]C[8*]", "[*]C[*]", "[3*]O[3*]"]


def test_precompute_properties():
    assert MOL_PROPERTIES_PATH.exists()


def test_filter_compound_init():
    compound_filter = CompoundFilter(MOL_PROPERTIES_PATH, fpSize=2048)
    assert compound_filter.len_BBs == len(n1_stock_subset)


@pytest.mark.parametrize(
    "fragment_smiles",
    [
        "[5*]N1CCC[C@@]1([13*])C",
        "[4*]CCN[5*]",
        "[4*]C[8*]",
        "[*]C[*]",
        "[3*]O[3*]",
    ],
)
def test_filter_compound(fragment_smiles):
    """The result from direct substructure match should be the same as after compound filtering."""
    compound_filter = CompoundFilter(MOL_PROPERTIES_PATH, fpSize=2048)
    substructure_matcher = SubstructureMatcher(set(n1_stock_subset))
    # without filter
    direct_valid_BBs = substructure_matcher.get_substructure_BBs(fragment_smiles)

    # with filter
    _, filtered_BBs = compound_filter.get_filtered_BBs(fragment_smiles)
    filtered_substructure_matcher = SubstructureMatcher(filtered_BBs, useChirality=True)
    filtered_valid_BBs = filtered_substructure_matcher.get_substructure_BBs(fragment_smiles)
    assert (
        direct_valid_BBs == filtered_valid_BBs
    ), f"Direct and filtered results differ for fragment SMILES: {fragment_smiles}"
    # assert len(direct_valid_BBs) == len(filtered_valid_BBs), f"Direct and filtered results differ for fragment SMILES: {fragment_smiles}"


if __name__ == "__main__":
    pytest.main()

    # import os
    # os.remove(PRECOMPUTE_PATH / "n1_stock_properties_first_500.json")
