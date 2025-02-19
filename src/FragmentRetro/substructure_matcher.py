import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Optional, cast

from rdkit import Chem

from FragmentRetro.utils.logging_config import logger
from FragmentRetro.utils.type_definitions import BBsType


class SubstructureMatcher:
    def __init__(
        self,
        BBs: BBsType = set(),
        useChirality: bool = True,
        parallelize: bool = False,
        num_cores: Optional[int] = None,
    ):
        """
        Initialize with a set of building blocks (BBs).

        Args:
            BBs: Set of building block SMILES strings.
            useChirality: Whether to match chirality.
            parallelize: Whether to enable parallel processing.
            num_cores: Number of CPU cores to use for parallel processing. If None, uses all available cores.
        """
        self.BBs = BBs
        self.useChirality = useChirality
        self.parallelize = parallelize
        # Use all available cores if not specified
        self.num_cores = num_cores if num_cores is not None else cpu_count()

    @staticmethod
    def convert_to_smarts(fragment_smiles: str) -> str:
        """
        Convert a fragment SMILES string into a SMARTS string with hydrogen counts
        explicitly specified, and convert dummy atoms into wildcards.

        Args:
            fragment_smiles: SMILES string of the fragment.

        Returns:
            SMARTS string with explicit hydrogen counts, handling dummy atoms.
        """
        mol = Chem.MolFromSmiles(fragment_smiles)
        # Add indices to atoms so that the `replace` function would
        # not replace atoms with the same atomic number
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {fragment_smiles}")
        smarts = Chem.MolToSmarts(mol)
        # Replace only specific indices
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                continue
            num_hydrogens = atom.GetTotalNumHs()
            idx = atom.GetAtomMapNum()
            smarts = smarts.replace(
                f"[#{atom.GetAtomicNum()}:{idx}]",
                f"[#{atom.GetAtomicNum()}H{num_hydrogens}:{idx}]",
            )
        # Remove indices
        smarts = re.sub(r":\d+\]", "]", smarts)
        # sub [#0] with *, which is a wildcard for any atom
        smarts = re.sub(r"\[\d*\#0\]", "*", smarts)
        return cast(str, smarts)

    @staticmethod
    def addH_to_wildcard_neighbors(fragment_smarts: str) -> str:
        """
        Adjust SMARTS by adding hydrogens to neighbors of wildcard atoms (*).

        Args:
            fragment_smarts: Fragment SMARTS string.

        Returns:
            Adjusted SMARTS string with explicit hydrogen counts for neighbors of wildcards.
        """
        # TODO: make sure this works properly and see if we want to merge this with convert_to_smarts

        # Convert SMARTS to molecule
        fragment_mol = Chem.MolFromSmarts(fragment_smarts)
        if fragment_mol is None:
            raise ValueError(f"Invalid SMARTS string: {fragment_smarts}")

        # Assign map indices to each atom
        for i, atom in enumerate(fragment_mol.GetAtoms()):
            atom.SetAtomMapNum(i + 1)

        # Convert molecule back to SMARTS with indices
        smarts_with_indices = Chem.MolToSmarts(fragment_mol)

        # Sanitize mol to avoid aromatic valency problems
        Chem.SanitizeMol(fragment_mol)

        for atom in fragment_mol.GetAtoms():
            addH = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "*":  # Neighbor is also wildcard
                    addH = True
                    break

            if addH:
                num_hydrogens = atom.GetTotalNumHs()
                idx = atom.GetAtomMapNum()

                # Update SMARTS with explicit hydrogen count
                smarts_with_indices = smarts_with_indices.replace(
                    # The '&' here is from `Chem.MolToSmarts(fragment_mol)`
                    f"[#{atom.GetAtomicNum()}&H{num_hydrogens}:{idx}]",
                    f"[#{atom.GetAtomicNum()}&H{num_hydrogens},#{atom.GetAtomicNum()}&H{num_hydrogens+1}]",
                )

        # Remove atom map indices
        adjusted_smarts = re.sub(r":\d+\]", "]", smarts_with_indices)
        return adjusted_smarts

    @staticmethod
    def is_strict_substructure(fragment_smiles: str, molecule_smiles: str, useChirality: bool = True) -> bool:
        """
        Check if the fragment is a strict substructure of the molecule. No extra atoms or
        branchings should be in the molecule other than the ones explicitly defined in the
        fragment, except for dummy atoms.

        Args:
            fragment_smiles: SMILES string of the fragment.
            molecule_smiles: SMILES string of the molecule.
            useChirality: whether to match chirality

        Returns:
            True if the fragment is a strict substructure of the molecule, False otherwise.
        """
        # Convert fragment SMILES to SMARTS pattern
        # hydrogen counts are explicitly specified to strictly match the fragment
        fragment_smarts = SubstructureMatcher.convert_to_smarts(fragment_smiles)
        # Add hydrogen atoms to wildcard neighbors to ensure hydrogen atoms can match with wildcard atoms
        fragment_smarts_withH = SubstructureMatcher.addH_to_wildcard_neighbors(fragment_smarts)

        # Convert molecule SMILES to RDKit molecule object
        # fragment_mol = Chem.MolFromSmarts(fragment_smarts)
        fragment_mol_withH = Chem.MolFromSmarts(fragment_smarts_withH)
        molecule_mol = Chem.MolFromSmiles(molecule_smiles)
        # to make sure hydrogen atoms can match with wildcard atoms
        molecule_mol = Chem.AddHs(molecule_mol)

        if molecule_mol is None:
            raise ValueError(f"Invalid SMILES string: {molecule_smiles}")
        return cast(bool, molecule_mol.HasSubstructMatch(fragment_mol_withH, useChirality=useChirality))

    def get_substructure_BBs(self, fragment: str) -> BBsType:
        """
        Get the set of building blocks (BBs) that the fragment matches.

        Args:
            fragment: SMILES string of the fragment.

        Returns:
            Set of building block SMILES strings that the fragment matches.
        """
        logger.info(f"[SubstructureMatcher] Matching fragment {fragment} to building blocks")
        if self.parallelize:
            logger.info(f"[SubstructureMatcher] Using {self.num_cores} cores for parallel processing")
            with ProcessPoolExecutor(max_workers=self.num_cores) as executor:
                future_to_bb = {
                    executor.submit(self.is_strict_substructure, fragment, bb, self.useChirality): bb for bb in self.BBs
                }

                strict_substructure_BBs = set()
                for future in as_completed(future_to_bb):
                    bb = future_to_bb[future]
                    try:
                        if future.result():  # If the result is True, add the building block to the set
                            strict_substructure_BBs.add(bb)
                    except Exception as exc:
                        logger.error(f"[SubstructureMatcher] Building block {bb} generated an exception: {exc}")
                        raise Exception("Parallel processing failed")
        else:
            # Fallback to single-threaded execution
            strict_substructure_BBs = set(
                bb for bb in self.BBs if self.is_strict_substructure(fragment, bb, self.useChirality)
            )
        logger.info(f"[SubstructureMatcher] Found {len(strict_substructure_BBs)} matching building")
        return strict_substructure_BBs
