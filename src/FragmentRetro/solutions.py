import itertools
from itertools import chain

from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.type_definitions import CombType, SolutionType


class RetrosynthesisSolution:
    def __init__(self, retrosynthesis: Retrosynthesis):
        self.retrosynthesis = retrosynthesis
        self.solutions: list[SolutionType] = []
        self.valid_combinations = list(chain.from_iterable(retrosynthesis.valid_combinations_dict.values()))
        self.num_fragments = retrosynthesis.fragmenter.num_fragments

    @staticmethod
    def get_solutions(valid_combinations: list[CombType], num_fragments: int) -> list[SolutionType]:
        """
        Generates all possible retrosynthesis solutions from a list of valid fragment combinations.

        A valid solution is defined as a combination of fragment lists that, when combined,
        include all the original fragments exactly once.

        Args:
            valid_combinations: A list of valid fragment combinations, where each combination
                is a list of fragment indices.
            num_fragments: The total number of fragments in the retrosynthesis.

        Returns:
            A list of retrosynthesis solutions. Each solution is a list of fragment lists.
            Returns an empty list if no solutions are found.
        """
        all_solutions: list[SolutionType] = []
        for r in range(1, num_fragments + 1):
            for solution in itertools.combinations(valid_combinations, r):
                flat_solution = [item for sublist in solution for item in sublist]
                if len(flat_solution) == num_fragments and set(flat_solution) == set(range(num_fragments)):
                    sorted_solution = sorted(solution, key=lambda s: sorted(s))
                    if sorted_solution not in all_solutions:
                        all_solutions.append(sorted_solution)

        return all_solutions

    def fill_solutions(self) -> None:
        """Fill the solutions list with all possible solutions."""
        self.solutions = self.get_solutions(self.valid_combinations, self.num_fragments)

    # def sort_solutions_by_len(self):
    #     return sorted(self.solutions, key=lambda x: len(x))

    # def visualize_solutions(self) -> None:
    #     pass
