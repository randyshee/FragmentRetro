from itertools import chain
from typing import Optional

from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.utils.logging_config import logger
from FragmentRetro.utils.type_definitions import CombType, SolutionType
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


class RetrosynthesisSolution:
    def __init__(self, retrosynthesis: Retrosynthesis):
        self.retrosynthesis = retrosynthesis
        self.fragmenter = retrosynthesis.fragmenter
        self.solutions: list[SolutionType] = []
        self.valid_combinations = list(chain.from_iterable(retrosynthesis.valid_combinations_dict.values()))
        self.valid_combinations = sorted(self.valid_combinations, key=len, reverse=True)
        self.num_fragments = retrosynthesis.fragmenter.num_fragments

    @staticmethod
    def _find_complementary_combinations(
        current_solution: list[CombType],
        remaining_fragments: set[int],
        start_index: int,
        valid_combinations: list[CombType],
        all_solutions: list[SolutionType],
    ) -> None:
        """
        Recursively finds complementary combinations to form a complete solution.

        Args:
            current_solution: The current partial solution.
            remaining_fragments: The set of fragments that still need to be covered.
            start_index: The index to start searching for combinations from.
            valid_combinations: The list of valid fragment combinations.
            all_solutions: The list to store complete solutions.
        """
        if not remaining_fragments:
            sorted_solution = sorted(current_solution, key=lambda s: sorted(s))
            if sorted_solution not in all_solutions:
                all_solutions.append(sorted_solution)
            return

        for j in range(start_index, len(valid_combinations)):
            comb2 = valid_combinations[j]
            if set(comb2).issubset(remaining_fragments):
                RetrosynthesisSolution._find_complementary_combinations(
                    current_solution + [comb2],
                    remaining_fragments - set(comb2),
                    j + 1,
                    valid_combinations,
                    all_solutions,
                )

    @staticmethod
    def get_solutions(
        valid_combinations: list[CombType], num_fragments: int, solution_cap: Optional[int] = None
    ) -> list[SolutionType]:
        """
        Generates all possible retrosynthesis solutions from a list of valid fragment combinations.

        A valid solution is defined as a combination of fragment lists that, when combined,
        include all the original fragments exactly once.

        Args:
            valid_combinations: A list of valid fragment combinations, where each combination
                is a tuple of fragment indices.
            num_fragments: The total number of fragments in the retrosynthesis.
            solution_cap: An optional integer specifying the maximum number of solutions to return.

        Returns:
            A list of retrosynthesis solutions. Each solution is a list of fragment lists.
            Returns an empty list if no solutions are found.
        """
        all_solutions: list[SolutionType] = []
        full_fragment_set = set(range(num_fragments))

        for i, comb1 in enumerate(valid_combinations):
            solution = [comb1]
            remaining_fragments = full_fragment_set - set(comb1)

            if not remaining_fragments:
                sorted_solution = sorted(solution, key=lambda s: sorted(s))
                if sorted_solution not in all_solutions:
                    all_solutions.append(sorted_solution)
                continue

            RetrosynthesisSolution._find_complementary_combinations(
                solution, remaining_fragments, i + 1, valid_combinations, all_solutions
            )
            if solution_cap and len(all_solutions) >= solution_cap:
                logger.info(f"[RetrosynthesisSolution] Solution count capped at {solution_cap}")
                return all_solutions[:solution_cap]

        return all_solutions

    def fill_solutions(self, solution_cap: Optional[int] = None) -> None:
        """Fill the solutions list with all possible solutions."""
        self.solutions = self.get_solutions(self.valid_combinations, self.num_fragments, solution_cap)

    def get_solution_smiles(self, solution: SolutionType) -> list[str]:
        """
        Retrieves the SMILES strings for each combination in a given solution.

        Args:
            solution: A retrosynthesis solution, which is a list of fragment combinations.
                      Each combination is a tuple of fragment indices.

        Returns:
            A list of SMILES strings, where each string corresponds to a fragment
            combination in the solution.
        """
        all_smiles = []
        for comb in solution:
            smiles = self.fragmenter.get_combination_smiles(comb)
            all_smiles.append(smiles)
        return all_smiles

    def visualize_solutions(
        self, solutions: list[SolutionType], molsPerRow: int = 3, subImgSize: tuple[float, float] = (200, 200)
    ) -> list[Image.Image]:
        """Visualizes a list of retrosynthesis solutions.

        Generates a list of images, where each image visualizes a single retrosynthesis
        solution. Each solution is displayed as a grid of molecules, with the SMILES
        strings of the molecules as legends.

        Args:
            solutions: A list of retrosynthesis solutions. Each solution is a list of
                fragment combinations (lists of fragment indices).
            molsPerRow: The number of molecules to display in each row of the grid.
            subImgSize: The size (width, height) of each molecule image in the grid.

        Returns:
            A list of PIL Image objects, each visualizing a retrosynthesis solution.
        """
        all_img = []
        for solution in solutions:
            logger.info(f"[RetrosynthesisSolution] Solution: {solution}")
            all_smiles = self.get_solution_smiles(solution)
            logger.info(f"[RetrosynthesisSolution] SMILES: {all_smiles}")
            # Convert SMILES to RDKit molecules
            mols = [Chem.MolFromSmiles(smiles) for smiles in all_smiles]
            # Draw molecules in a grid
            legends = [f"{comb}: {smiles}" for comb, smiles in zip(solution, all_smiles, strict=False)]
            img = Draw.MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends)
            all_img.append(img)
        return all_img
