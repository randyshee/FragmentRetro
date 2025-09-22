import re
from multiprocessing import Pool, cpu_count
from typing import cast

from rdkit import Chem
from rdkit.Chem.rdchem import Atom

from fragmentretro.typing import BBsType
from fragmentretro.utils.logging_config import logger


class SubstructureMatcher:
    def __init__(
        self,
        BBs: BBsType,
        useChirality: bool = True,
        parallelize: bool = False,
        num_cores: int | None = None,
        core_factor: int = 10,
    ):
        """
        Initialize with a set of building blocks (BBs).

        Args:
            BBs: Set of building block SMILES strings.
            useChirality: Whether to match chirality.
            parallelize: Whether to enable parallel processing.
            num_cores: Number of CPU cores to use for parallel processing. If None, uses all available cores.
            core_factor: Factor to determine when to use parallel processing. If the number of BBs is greater than
                         `num_cores * core_factor`, parallel processing will be used.
        """
        self.BBs = BBs
        self.useChirality = useChirality
        self.parallelize = parallelize
        # Use all available cores if not specified
        self.num_cores = min(num_cores if num_cores is not None else cpu_count(), cpu_count())
        self.core_factor = core_factor

    @staticmethod
    def has_dummy_neighbor(atom: Atom) -> bool:
        dummy_neighbor = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == "*":  # Neighbor is wildcard
                dummy_neighbor = True
                break
        return dummy_neighbor

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
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 0:
                continue
            num_hydrogens = atom.GetTotalNumHs()
            idx = atom.GetAtomMapNum()
            if num_hydrogens == 0:
                smarts = re.sub(rf"\[\#{atomic_num}(@*):{idx}\]", rf"[#{atomic_num}\1&H{num_hydrogens}:{idx}]", smarts)
            else:
                smarts = re.sub(rf"\[\#{atomic_num}:{idx}\]", rf"[#{atomic_num}&H{num_hydrogens}:{idx}]", smarts)
        # Remove indices
        smarts = re.sub(r":\d+\]", "]", smarts)
        # sub [#0] with *, which is a wildcard for any atom
        smarts = re.sub(r"\[\d*\#0\]", "*", smarts)
        return smarts

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
            dummy_neighbor = SubstructureMatcher.has_dummy_neighbor(atom)

            if dummy_neighbor:
                atomic_num = atom.GetAtomicNum()
                num_hydrogens = atom.GetTotalNumHs()
                idx = atom.GetAtomMapNum()

                logger.debug("smarts_with_indices: %s", smarts_with_indices)
                logger.debug("updating atom idx: %s", idx)
                # Update SMARTS with explicit hydrogen count
                # remove chirality information (in re group \1) for atoms with dummy neighbors
                smarts_with_indices = re.sub(
                    rf"\[\#{atomic_num}(@*)&H{num_hydrogens}:{idx}\]",
                    rf"[#{atomic_num}&H{num_hydrogens},#{atomic_num}&H{num_hydrogens + 1}:{idx}]",
                    smarts_with_indices,
                )
                logger.debug("smarts_with_indices: %s", smarts_with_indices)

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
        logger.debug(f"[SubstructureMatcher] Matching fragment {fragment} to building blocks")
        if self.parallelize and len(self.BBs) >= self.num_cores * self.core_factor:
            logger.debug(f"[SubstructureMatcher] Using {self.num_cores} cores for parallel processing")
            with Pool(processes=self.num_cores) as pool:
                results = pool.starmap(
                    self.is_strict_substructure, [(fragment, bb, self.useChirality) for bb in self.BBs]
                )
                strict_substructure_BBs = set(bb for bb, result in zip(self.BBs, results, strict=False) if result)

        else:
            # Fallback to single-threaded execution
            strict_substructure_BBs = set(
                bb for bb in self.BBs if self.is_strict_substructure(fragment, bb, self.useChirality)
            )
        logger.debug(f"[SubstructureMatcher] Found {len(strict_substructure_BBs)} matching building")
        return strict_substructure_BBs
