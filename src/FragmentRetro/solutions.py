import itertools
from itertools import chain

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

from FragmentRetro.logging_config import logger
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.type_definitions import CombType, SolutionType


class RetrosynthesisSolution:
    def __init__(self, retrosynthesis: Retrosynthesis):
        self.retrosynthesis = retrosynthesis
        self.fragmenter = retrosynthesis.fragmenter
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
            logger.info(f"Solution: {solution}")
            all_smiles = []
            for comb in solution:
                smiles = self.fragmenter._get_combination_smiles(comb)
                all_smiles.append(smiles)
            logger.info(f"SMILES: {all_smiles}")
            # Convert SMILES to RDKit molecules
            mols = [Chem.MolFromSmiles(smiles) for smiles in all_smiles]
            # Draw molecules in a grid
            img = Draw.MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=all_smiles)
            all_img.append(img)
        return all_img
