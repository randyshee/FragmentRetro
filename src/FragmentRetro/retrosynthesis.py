from pathlib import Path
from typing import Optional, cast

from FragmentRetro.fragmenter_base import Fragmenter
from FragmentRetro.substructure_matcher import SubstructureMatcher
from FragmentRetro.utils.filter_compound import CompoundFilter
from FragmentRetro.utils.helpers import replace_dummy_atoms_regex
from FragmentRetro.utils.logging_config import logger
from FragmentRetro.utils.type_definitions import (
    BBsType,
    CombBBsDictType,
    CombType,
    StageCombDictType,
)


class Retrosynthesis:
    def __init__(
        self, fragmenter: Fragmenter, original_BBs: Optional[BBsType] = None, mol_properties_path: Optional[Path] = None
    ):
        self.fragmenter = fragmenter
        self.num_fragments = fragmenter.num_fragments
        self.valid_combinations_dict: StageCombDictType = {}  # store valid combs for each stage
        self.invalid_combinations_dict: StageCombDictType = {}  # store invalid combs for each stage
        self.comb_bbs_dict: CombBBsDictType = {}  # store valid BBs for fragment combs

        if original_BBs is not None and mol_properties_path is not None:
            logger.warn("Both original_BBs and mol_properties_path are provided. " "Will be using mol_properties_path.")
        elif original_BBs is None and mol_properties_path is None:
            logger.critical("Either original_BBs or mol_properties_path must be provided.")
        if mol_properties_path is not None:
            self.use_filter = True
            self.compound_filter = CompoundFilter(mol_properties_path)
        else:
            self.use_filter = False
            self.original_BBs = original_BBs

    def _check_effective_comb(self, comb: CombType) -> bool:
        """Check if a combination is effective.

        A combination is considered effective if it does not contain any
        invalid combinations from the previous stage. Note that only effective
        combinations could become a valid combination.

        Args:
            comb: The combination to check, represented as a tuple of fragment indices.

        Returns:
            True if the combination is effective, False otherwise.
        """
        len_comb = len(comb)
        len_minus_one_invalid_combs = self.invalid_combinations_dict.get(len_comb - 1, [])
        for invalid_comb in len_minus_one_invalid_combs:
            if set(invalid_comb).issubset(set(comb)):
                return False
        return True

    def _get_possible_BBs_for_comb(self, comb: CombType) -> BBsType:
        """Get possible building blocks for a given combination of fragments.

        For a combination of length 1, the original building blocks are returned.
        For combinations of length greater than 1, the building blocks are
        retrieved by intersecting the building blocks of valid combinations from
        the previous stage that are subsets of the current combination.

        Args:
            comb: The combination to get building blocks for, represented as a
                tuple of fragment indices.

        Returns:
            A set of SMILES strings representing the building blocks for the
            given combination.
        """
        len_comb = len(comb)
        if len_comb == 1:
            if self.use_filter:
                comb_smiles = self.fragmenter.get_combination_smiles(comb)
                no_dummy_smiles = replace_dummy_atoms_regex(comb_smiles)
                return self.compound_filter.get_filtered_BBs(no_dummy_smiles)
            else:
                return cast(BBsType, self.original_BBs)
        len_minus_one_valid_combs = self.valid_combinations_dict.get(len_comb - 1, [])
        BBs: BBsType = set()
        for valid_comb in len_minus_one_valid_combs:
            if set(valid_comb).issubset(set(comb)):
                if not BBs:
                    BBs = self.comb_bbs_dict[valid_comb]
                else:
                    BBs = BBs.intersection(self.comb_bbs_dict[valid_comb])
        return BBs

    def _retro_stage(self, stage: int) -> None:
        """Perform retrosynthesis for a single stage.

        This method performs the retrosynthesis process for a given stage. It
        generates combinations of fragments, checks their effectiveness,
        identifies valid combinations based on building block matching, and
        stores the results.

        Args:
            stage: The current retrosynthesis stage (an integer).
        """
        self.valid_combinations_dict[stage] = []
        # get fragment comb for stage
        combs = list(self.fragmenter.get_length_n_combinations(stage))
        logger.info(f"Stage {stage}: {len(combs)} combinations")
        # check invalid comb and filter out effective comb
        effective_combs, invalid_combs = [], []
        for comb in combs:
            if self._check_effective_comb(comb):
                effective_combs.append(comb)
            else:
                invalid_combs.append(comb)
        self.invalid_combinations_dict[stage] = invalid_combs
        logger.info(f"Stage {stage}: {len(effective_combs)} effective combinations")

        for comb in effective_combs:
            fragment_smiles = self.fragmenter.get_combination_smiles(comb)
            # get building blocks for comb
            possible_comb_BBs = self._get_possible_BBs_for_comb(comb)
            comb_matcher = SubstructureMatcher(possible_comb_BBs)
            valid_BBs = comb_matcher.get_substructure_BBs(fragment_smiles)
            # store valid comb and BBs
            if len(valid_BBs) > 0:
                self.valid_combinations_dict[stage].append(comb)
                self.comb_bbs_dict[comb] = valid_BBs
            else:
                self.invalid_combinations_dict[stage].append(comb)
        if stage == 1 and self.use_filter:
            del self.compound_filter  # save memory
        logger.info(f"Stage {stage}: {len(self.valid_combinations_dict[stage])} valid combinations")
        logger.info(f"Stage {stage}: {len(self.invalid_combinations_dict[stage])} invalid combinations")

    def fragment_retrosynthesis(self) -> StageCombDictType:
        """Perform retrosynthesis on the molecule.

        This method orchestrates the retrosynthesis process by iteratively
        applying the `_retro_stage` method for each stage, starting from stage 1
        up to the total number of fragments.

        Returns:
            A dictionary containing the valid combinations for each stage.
            The keys are the stage numbers (integers), and the values are lists
            of valid combinations (tuples of fragment indices).
        """
        for stage in range(1, self.num_fragments + 1):
            self._retro_stage(stage)
        return self.valid_combinations_dict
