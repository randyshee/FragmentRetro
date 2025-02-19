"""Tests for substructure matcher."""

from pathlib import Path

import pytest

from FragmentRetro.substructure_matcher import SubstructureMatcher

DATA_PATH = Path(__file__).parent.parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"

TEST_CASES_FOR_CONVERT_TO_SMARTS = [
    {
        "case_number": 1,
        "fragment_smiles": "[1*]C([6*])=O",
        "expected_smarts": "*-[#6H0](-*)=[#8H0]",
        "description": "Without hydrogen and two dummy atoms",
    },
    {
        "case_number": 2,
        "fragment_smiles": "[14*]c1nc2c(C)c(Cl)ccc2[nH]1",
        "expected_smarts": "*-[#6H0]1:[#7H0]:[#6H0]2:[#6H0](-[#6H3]):[#6H0](-[#17H0]):[#6H1]:[#6H1]:[#6H0]:2:[#7H]:1",
        "description": "Two rings",
    },
    {
        "case_number": 3,
        "fragment_smiles": "[16*]c1ccc([16*])c([16*])c1",
        "expected_smarts": "*-[#6H0]1:[#6H1]:[#6H1]:[#6H0](-*):[#6H0](-*):[#6H1]:1",
        "description": "Three dummy atoms and double-digit indices",
    },
    {
        "case_number": 4,
        "fragment_smiles": "[3*]OC",
        "expected_smarts": "*-[#8H0]-[#6H3]",
        "description": "Three hydrogen atoms on carbon",
    },
    {
        "case_number": 5,
        "fragment_smiles": "[5*]N1CCC[C@]1([13*])C",
        "expected_smarts": "*-[#7H0]1-[#6H2]-[#6H2]-[#6H2]-[#6@]-1(-*)-[#6H3]",
        "description": "Chiral",
    },
    {
        "case_number": 6,
        "fragment_smiles": "[9*]n1nccn1",
        "expected_smarts": "*-[#7H0]1:[#7H0]:[#6H1]:[#6H1]:[#7H0]:1",
        "description": "Aromatic",
    },
    {
        "case_number": 7,
        "fragment_smiles": "[*]CC",
        "expected_smarts": "*-[#6H2]-[#6H3]",
        "description": "Dummy atoms without isotopic specification",
    },
    {
        "case_number": 8,
        "fragment_smiles": "*CC",
        "expected_smarts": "*-[#6H2]-[#6H3]",
        "description": "Dummy atoms without isotopic specification and brackets",
    },
    {
        "case_number": 9,
        "fragment_smiles": "CCNC",
        "expected_smarts": "[#6H3]-[#6H2]-[#7H1]-[#6H3]",
        "description": "SMILES without dummy atoms",
    },
    # TODO: Add test cases for molecules with charges
]

TEST_CASES_FOR_IS_STRICT_SUBSTRUCTURE = [
    {
        "case_number": 1,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles": "CCCN",
        "expected": True,
        "description": "One atom branching out from dummy atom",
    },
    {
        "case_number": 2,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles": "CCCCN",
        "expected": True,
        "description": "Two atoms branching out from dummy atom",
    },
    {
        "case_number": 3,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles": "CCCNC",
        "expected": False,
        "description": "One atom branching out from non-dummy atom",
    },
    {
        "case_number": 4,
        "fragment_smiles": "[4*]CCN[5*]",
        "molecule_smiles": "CCCNC",
        "expected": True,
        "description": "Each atom branching out from the two non-dummy atoms",
    },
    {
        "case_number": 5,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles": "NCCCN",
        "expected": True,
        "description": "Non-carbon atom branching out from dummy atom",
    },
    {
        "case_number": 6,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles": "CCN",
        "expected": True,
        "description": "No atom branching out from dummy atom should also be true",
    },
    {
        "case_number": 7,
        "fragment_smiles": "[4*]CCN[5*]",
        "molecule_smiles": "CCN",
        "expected": True,
        "description": "No atom branching out from two dummy atoms should also be true",
    },
    {
        "case_number": 8,
        "fragment_smiles": "[4*]CCN[5*]",
        "molecule_smiles": "CCNC",
        "expected": True,
        "description": "No atom branching out from just one dummy atom should also be true",
    },
    {
        "case_number": 8,
        "fragment_smiles": "[9*]n1cnc2ncc([14*])nc21",
        "molecule_smiles": "Cc1ccc(S(=O)(=O)n2ccc3nc(NN)cnc32)cc1",
        "expected": False,
        "description": "Wrong atom at the neighbor of a dummy atom should be false",
    },
    {
        "case_number": 10,
        "fragment_smiles": "[16*]c1ccc([16*])c([16*])c1",
        "molecule_smiles": "Brc1ccc(Br)cc1",
        "expected": True,
        "description": "No atom branching out from just two dummy atom should also be true",
    },
    {
        "case_number": 11,
        "fragment_smiles": "CCNC",
        "molecule_smiles": "CCNC",
        "expected": True,
        "description": "SMILES without * should match normally",
    },
    {
        "case_number": 12,
        "fragment_smiles": "CCNC",
        "molecule_smiles": "CCCNC",
        "expected": False,
        "description": "SMILES without * should not match with other SMILES",
    },
    {
        "case_number": 13,
        "fragment_smiles": "[5*]N1CCC[C@@]1([13*])C",
        "molecule_smiles": "BrN1CCC[C@@]1(Br)C",
        "expected": True,
        "description": "Chiral case True",
    },
    {
        "case_number": 14,
        "fragment_smiles": "[5*]N1CCC[C@@]1([13*])C",
        "molecule_smiles": "BrN1CCC[C@]1(Br)C",
        "expected": False,
        "description": "Chiral case False",
    },
    # TODO: Add test cases for molecules with charges
]


@pytest.mark.parametrize(
    "case_number, fragment_smiles, expected_smarts, description",
    [
        (
            tc["case_number"],
            tc["fragment_smiles"],
            tc["expected_smarts"],
            tc["description"],
        )
        for tc in TEST_CASES_FOR_CONVERT_TO_SMARTS
    ],
)
def test_convert_to_smarts(case_number, fragment_smiles, expected_smarts, description):
    result = SubstructureMatcher.convert_to_smarts(fragment_smiles)
    assert result == expected_smarts, f"Case {case_number} failed: {description}. Fragment SMILES: {fragment_smiles}"


@pytest.mark.parametrize(
    "case_number, fragment_smiles, molecule_smiles, expected, description",
    [
        (
            tc["case_number"],
            tc["fragment_smiles"],
            tc["molecule_smiles"],
            tc["expected"],
            tc["description"],
        )
        for tc in TEST_CASES_FOR_IS_STRICT_SUBSTRUCTURE
    ],
)
def test_is_strict_substructure(case_number, fragment_smiles, molecule_smiles, expected, description):
    result = SubstructureMatcher.is_strict_substructure(fragment_smiles, molecule_smiles)
    assert result == expected, f"Case {case_number} failed: {description}"


@pytest.mark.parametrize(
    "fragment_smiles",
    [
        "[5*]N1CCC[C@@]1([13*])C",
        "[4*]CCN[5*]",
        "[4*]C[8*]",
        "[*]C[*]",
        "[3*]O[3*]",
        "[*]c1ccc(Nc2ncn(-c3ccnc(N4CC(C)N(C(C)=O)C(C)C4)c3)n2)cc1",
    ],
)
def test_parallel_get_substructure_BBs(fragment_smiles):
    with open(PAROUTES_PATH / "n1-stock.txt", "r") as f:
        n1_stock = [line.strip() for line in f.readlines()]
    n1_stock_subset = set(n1_stock[:500])
    no_parallel_matcher = SubstructureMatcher(n1_stock_subset, parallelize=False)
    parallel_matcher = SubstructureMatcher(n1_stock_subset, parallelize=True, num_cores=5, core_factor=5)
    assert no_parallel_matcher.get_substructure_BBs(fragment_smiles) == parallel_matcher.get_substructure_BBs(
        fragment_smiles
    ), f"Parallel and non-parallel results differ for fragment SMILES: {fragment_smiles}"


if __name__ == "__main__":
    pytest.main()
