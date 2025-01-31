from FragmentRetro.fragmenter_base import Fragmenter
from FragmentRetro.substructure_matcher import SubstructureMatcher
from FragmentRetro.type_definitions import CombBBsDictType, StageCombDictType


class Retrosynthesis:
    def __init__(self, fragmenter: Fragmenter, substructure_matcher: SubstructureMatcher):
        self.fragmenter = fragmenter
        self.matcher = substructure_matcher
        self.valid_combinations: StageCombDictType = {}
        self.invalid_combinations: StageCombDictType = {}
        self.comb_bbs: CombBBsDictType = {}

    # def fragment_retrosynthesis(self, stage: int = 1) -> list[set[list[int]]]:
    #     """
    #     Perform retrosynthesis analysis.

    #     Args:
    #         stage: The retrosynthesis stage.

    #     Returns:
    #         List of sets of fragment combinations.
    #     """
    #     if stage == 1:
    #         fragments = self.fragmenter._get_initial_fragments()
    #         for fragment in fragments:
    #             if not self.fragmenter.substructure_matcher.is_substructure_BBs(fragment, self.fragmenter.BBs):
    #                 logging.info(f"Fragment {fragment} is not a building block.")
    #                 return []
    #     n_combinations = self.fragmenter._get_length_n_combinations(n=stage)
    #     combination_smiles = [self.fragmenter._get_combination_smiles(comb) for comb in n_combinations]

    #     # have to record valid and invalid combination (comb of node IDs) not the smiles
    #     pass
