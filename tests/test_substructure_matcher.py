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
        "expected_smarts": "*-[#6&H0](-*)=[#8&H0]",
        "description": "Without hydrogen and two dummy atoms",
    },
    {
        "case_number": 2,
        "fragment_smiles": "[14*]c1nc2c(C)c(Cl)ccc2[nH]1",
        "expected_smarts": "*-[#6&H0]1:[#7&H0]:[#6&H0]2:[#6&H0](-[#6&H3]):[#6&H0](-[#17&H0]):[#6&H1]:[#6&H1]:[#6&H0]:2:[#7H]:1",
        "description": "Two rings",
    },
    {
        "case_number": 3,
        "fragment_smiles": "[16*]c1ccc([16*])c([16*])c1",
        "expected_smarts": "*-[#6&H0]1:[#6&H1]:[#6&H1]:[#6&H0](-*):[#6&H0](-*):[#6&H1]:1",
        "description": "Three dummy atoms and double-digit indices",
    },
    {
        "case_number": 4,
        "fragment_smiles": "[3*]OC",
        "expected_smarts": "*-[#8&H0]-[#6&H3]",
        "description": "Three hydrogen atoms on carbon",
    },
    {
        "case_number": 5,
        "fragment_smiles": "[5*]N1CCC[C@]1([13*])C",
        "expected_smarts": "*-[#7&H0]1-[#6&H2]-[#6&H2]-[#6&H2]-[#6@]-1(-*)-[#6&H3]",
        "description": "Chiral",
    },
    {
        "case_number": 6,
        "fragment_smiles": "[9*]n1nccn1",
        "expected_smarts": "*-[#7&H0]1:[#7&H0]:[#6&H1]:[#6&H1]:[#7&H0]:1",
        "description": "Aromatic",
    },
    {
        "case_number": 7,
        "fragment_smiles": "[*]CC",
        "expected_smarts": "*-[#6&H2]-[#6&H3]",
        "description": "Dummy atoms without isotopic specification",
    },
    {
        "case_number": 8,
        "fragment_smiles": "*CC",
        "expected_smarts": "*-[#6&H2]-[#6&H3]",
        "description": "Dummy atoms without isotopic specification and brackets",
    },
    {
        "case_number": 9,
        "fragment_smiles": "CCNC",
        "expected_smarts": "[#6&H3]-[#6&H2]-[#7&H1]-[#6&H3]",
        "description": "SMILES without dummy atoms",
    },
    {
        "case_number": 10,
        "fragment_smiles": "C[C@@H]([*])(N)",
        "expected_smarts": "[#6&H3]-[#6@@H](-*)-[#7&H2]",
        "description": "Chirality and hydrogen after atomic number",
    },
    # TODO: Add test cases for molecules with charges
]

TEST_CASES_ADDH_TO_WILDCARD_NEIGHBORS = [
    {
        "case_number": 1,
        "fragment_smarts": "*-[#7&H0]1-[#6&H2]-[#6&H2]-[#6&H2]-[#6@&H0]-1(-*)-[#6&H3]",
        "expected_smarts": "[*]-[#7&H0,#7&H1]1-[#6&H2]-[#6&H2]-[#6&H2]-[#6&H0,#6&H1]-1(-[*])-[#6&H3]",
        "description": "from fragment smiles [5*]N1CCC[C@]1([13*])C",
    },
    {
        "case_number": 2,
        "fragment_smarts": "[#6&H3]-[#6@@H&H1](-*)-[#7&H2]",
        "expected_smarts": "[#6&H3]-[#6H&H1,#6H&H2](-[*])-[#7&H2]",
        "description": "from fragment smiles C[C@@H]([*])(N)",
    },
]

TEST_CASES_FOR_IS_STRICT_SUBSTRUCTURE = [
    {
        "case_number": 1,
        "fragment_smiles": "[4*]CCN",
        "molecule_smiles_list": ["CCCN", "CCCCN", "CCCNC", "NCCCN", "CCN"],
        "expected_list": [True, True, False, True, True],
        "descriptions": [
            "One atom branching out from dummy atom",
            "Two atoms branching out from dummy atom",
            "One atom branching out from non-dummy atom",
            "Non-carbon atom branching out from dummy atom",
            "No atom branching out from dummy atom should also be true",
        ],
    },
    {
        "case_number": 2,
        "fragment_smiles": "[4*]CCN[5*]",
        "molecule_smiles_list": ["CCCNC", "CCN", "CCNC"],
        "expected_list": [True, True, True],
        "descriptions": [
            "Each atom branching out from the two non-dummy atoms",
            "No atom branching out from two dummy atoms should also be true",
            "No atom branching out from just one dummy atom should also be true",
        ],
    },
    {
        "case_number": 3,
        "fragment_smiles": "[9*]n1cnc2ncc([14*])nc21",
        "molecule_smiles_list": ["Cc1ccc(S(=O)(=O)n2ccc3nc(NN)cnc32)cc1"],
        "expected_list": [False],
        "descriptions": ["Wrong atom at the neighbor of a dummy atom should be false"],
    },
    {
        "case_number": 4,
        "fragment_smiles": "[16*]c1ccc([16*])c([16*])c1",
        "molecule_smiles_list": ["Brc1ccc(Br)cc1"],
        "expected_list": [True],
        "descriptions": ["No atom branching out from just two dummy atom should also be true"],
    },
    {
        "case_number": 5,
        "fragment_smiles": "CCNC",
        "molecule_smiles_list": ["CCNC", "CCCNC"],
        "expected_list": [True, False],
        "descriptions": [
            "SMILES without * should match normally",
            "SMILES without * should not match with other SMILES",
        ],
    },
    {
        "case_number": 6,
        "fragment_smiles": "[5*]N1CCC[C@@]1(Br)C",
        "molecule_smiles_list": ["BrN1CCC[C@@]1(Br)C", "BrN1CCC[C@]1(Br)C"],
        "expected_list": [True, False],
        "descriptions": [
            "Chiral case True",
            "Chiral case False",
        ],
    },
    {
        "case_number": 7,
        "fragment_smiles": "[5*]N1CCC[C@@]1([13*])C",
        "molecule_smiles_list": [
            "BrN1CCC[C@@]1(Br)C",
            "BrN1CCC[C@]1(Br)C",
            "BrN1CCC[C@@]1([H])C",
            "BrN1CCC[C@]1([H])C",
        ],
        "expected_list": [True, True, True, True],
        "descriptions": [
            "both chirality cases should be true",
            "both chirality cases should be true",
            "Hydrogen atom being the dummy atom at the neighbor of a chiral atom should be true",
            "both chirality cases should be true for hydrogen as the dummy atom",
        ],
    },
    {
        "case_number": 8,
        "fragment_smiles": "C[C@@H]([*])(N)",
        "molecule_smiles_list": ["CCN", "C[C@@H](Br)(N)", "C[C@H](Br)(N)", "C[C@@H](N)O", "C[C@H](N)O"],
        "expected_list": [True, True, True, True, True],
        "descriptions": [
            "Hydrogen atom being the dummy atom at the neighbor of a chiral atom should be true",
            "both chirality cases should be true",
            "both chirality cases should be true",
            "both chirality cases should be true for hydrogen as the dummy atom",
            "both chirality cases should be true for hydrogen as the dummy atom",
        ],
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
    "case_number, fragment_smarts, expected_smarts, description",
    [
        (
            tc["case_number"],
            tc["fragment_smarts"],
            tc["expected_smarts"],
            tc["description"],
        )
        for tc in TEST_CASES_ADDH_TO_WILDCARD_NEIGHBORS
    ],
)
def test_addH_to_wildcard_neighbors(case_number, fragment_smarts, expected_smarts, description):
    result = SubstructureMatcher.addH_to_wildcard_neighbors(fragment_smarts)
    assert result == expected_smarts, f"Case {case_number} failed: {description}"


@pytest.mark.parametrize(
    "case_number, fragment_smiles, molecule_smiles_list, expected_list, descriptions",
    [
        (
            tc["case_number"],
            tc["fragment_smiles"],
            tc["molecule_smiles_list"],
            tc["expected_list"],
            tc["descriptions"],
        )
        for tc in TEST_CASES_FOR_IS_STRICT_SUBSTRUCTURE
    ],
)
def test_is_strict_substructure(case_number, fragment_smiles, molecule_smiles_list, expected_list, descriptions):
    for i in range(len(molecule_smiles_list)):
        result = SubstructureMatcher.is_strict_substructure(fragment_smiles, molecule_smiles_list[i])
        assert result == expected_list[i], f"Case {case_number} failed: {descriptions[i]}"


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
