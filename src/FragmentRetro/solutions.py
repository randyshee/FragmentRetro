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
        # TODO: modify this method
        # Generate all combinations of valid combinations
        all_solutions: list[SolutionType] = []
        for r in range(1, num_fragments + 1):
            for combination in itertools.combinations(valid_combinations, r):
                # flatten the combination
                flat_combination = [item for sublist in combination for item in sublist]
                if set(flat_combination) == set(range(num_fragments)):
                    all_solutions.append(list(combination))

        return all_solutions

    def fill_solutions(self) -> None:
        """Fill the solutions list with all possible solutions."""
        self.solutions = self.get_solutions(self.valid_combinations, self.num_fragments)

    # def sort_solutions_by_len(self):
    #     return sorted(self.solutions, key=lambda x: len(x))

    # def visualize_solutions(self) -> None:
    #     pass
