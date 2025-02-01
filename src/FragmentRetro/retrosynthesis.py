from FragmentRetro.fragmenter_base import Fragmenter
from FragmentRetro.type_definitions import CombBBsDictType, StageCombDictType

# from FragmentRetro.substructure_matcher import SubstructureMatcher


class Retrosynthesis:
    def __init__(self, fragmenter: Fragmenter):
        self.fragmenter = fragmenter
        self.valid_combinations_dict: StageCombDictType = {}
        self.invalid_combinations_dict: StageCombDictType = {}
        self.comb_bbs_dict: CombBBsDictType = {}

    # def fragment_retrosynthesis(self, stage: int = 1) -> list[set[list[int]]]:
    #     """
    #     Perform retrosynthesis analysis.

    #     Args:
    #         stage: The retrosynthesis stage.

    #     Returns:
    #         List of sets of fragment combinations.
    #     """
    #     matcher = substructure_matcher # new matcher for each stage?
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
