"""Tests for solutions."""

import pytest

from FragmentRetro.solutions import RetrosynthesisSolution

TEST_CASES_FOR_GET_SOLUTIONS = [
    {
        "case_number": 1,
        "valid_combinations": [(0,), (1,), (2,), (0, 1), (1, 2), (0, 2), (0, 1, 2)],
        "num_fragments": 3,
        "expected_solution": [[(0, 1, 2)], [(0,), (1, 2)], [(0, 2), (1,)], [(0, 1), (2,)], [(0,), (1,), (2,)]],
        "description": "All possible combinations",
    },
    {
        "case_number": 2,
        "valid_combinations": [(0,), (1,), (2,), (0, 1), (1, 2)],
        "num_fragments": 3,
        "expected_solution": [[(0,), (1, 2)], [(0, 1), (2,)], [(0,), (1,), (2,)]],
        "description": "Some possible combinations",
    },
    {
        "case_number": 3,
        "valid_combinations": [(0,), (1,), (2,)],
        "num_fragments": 3,
        "expected_solution": [[(0,), (1,), (2,)]],
        "description": "Only one possible combination",
    },
    {
        "case_number": 4,
        "valid_combinations": [(0,), (1,), (2,), (3,), (0, 1), (1, 2), (0, 1, 2)],
        "num_fragments": 4,
        "expected_solution": [[(0, 1, 2), (3,)], [(0,), (1, 2), (3,)], [(0, 1), (2,), (3,)], [(0,), (1,), (2,), (3,)]],
        "description": "Some possible combinations for 4 fragments",
    },
    {
        "case_number": 5,
        "valid_combinations": [(0,), (1,), (2,), (3,), (4,), (0, 1), (1, 2), (0, 1, 2), (0, 1, 2, 3)],
        "num_fragments": 5,
        "expected_solution": [
            [(0, 1, 2, 3), (4,)],
            [(0, 1, 2), (3,), (4,)],
            [(0,), (1, 2), (3,), (4,)],
            [(0, 1), (2,), (3,), (4,)],
            [(0,), (1,), (2,), (3,), (4,)],
        ],
        "description": "Some possible combinations for 4 fragments",
    },
]


@pytest.mark.parametrize(
    "case_number, valid_combinations, num_fragments, expected_solution, description",
    [
        (
            tc["case_number"],
            tc["valid_combinations"],
            tc["num_fragments"],
            tc["expected_solution"],
            tc["description"],
        )
        for tc in TEST_CASES_FOR_GET_SOLUTIONS
    ],
)
def test_get_solutions(case_number, valid_combinations, num_fragments, expected_solution, description):
    result = RetrosynthesisSolution.get_solutions(valid_combinations, num_fragments)
    assert (
        len(result) == len(expected_solution)
    ), f"Case {case_number} failed: {description}. Length of results differs. Expected {len(expected_solution)}, got {len(result)}"
    for sol in expected_solution:
        assert sol in result, f"Case {case_number} failed: {description}. Expected solution {sol} not found in results."


if __name__ == "__main__":
    pytest.main()
