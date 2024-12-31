"""Tests for substructure matcher."""

import pytest

from FragmentRetro.substructure_matcher import SubstructureMatcher

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
    assert (
        result == expected_smarts
    ), f"Case {case_number} failed: {description}. Fragment SMILES: {fragment_smiles}"


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
def test_is_strict_substructure(
    case_number, fragment_smiles, molecule_smiles, expected, description
):
    result = SubstructureMatcher.is_strict_substructure(
        fragment_smiles, molecule_smiles
    )
    assert result == expected, f"Case {case_number} failed: {description}"


if __name__ == "__main__":
    pytest.main()
