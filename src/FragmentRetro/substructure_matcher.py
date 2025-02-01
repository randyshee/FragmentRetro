import re
from typing import cast

from rdkit import Chem


class SubstructureMatcher:
    def __init__(self, BBs: set[str] = set()):
        """
        Initialize with a set of building blocks (BBs).

        Args:
            BBs: Set of building block SMILES strings.
        """
        self.BBs = BBs

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
                    f"[#{atom.GetAtomicNum()}H{num_hydrogens+1}:{idx}]",
                )

        # Remove atom map indices
        adjusted_smarts = re.sub(r":\d+\]", "]", smarts_with_indices)
        return adjusted_smarts

    @staticmethod
    def is_strict_substructure(fragment_smiles: str, molecule_smiles: str) -> bool:
        """
        Check if the fragment is a strict substructure of the molecule. No extra atoms or
        branchings should be in the molecule other than the ones explicitly defined in the
        fragment, except for dummy atoms.

        Args:
            fragment_smiles: SMILES string of the fragment.
            molecule_smiles: SMILES string of the molecule.

        Returns:
            True if the fragment is a strict substructure of the molecule, False otherwise.
        """
        # Convert fragment SMILES to SMARTS pattern
        # hydrogen counts are explicitly specified to strictly match the fragment
        fragment_smarts = SubstructureMatcher.convert_to_smarts(fragment_smiles)
        # Add hydrogen atoms to wildcard neighbors to ensure hydrogen atoms can match with wildcard atoms
        fragment_smarts_withH = SubstructureMatcher.addH_to_wildcard_neighbors(fragment_smarts)

        # Convert molecule SMILES to RDKit molecule object
        fragment_mol = Chem.MolFromSmarts(fragment_smarts)
        fragment_mol_withH = Chem.MolFromSmarts(fragment_smarts_withH)
        molecule_mol = Chem.MolFromSmiles(molecule_smiles)
        # to make sure hydrogen atoms can match with wildcard atoms
        molecule_mol = Chem.AddHs(molecule_mol)

        if molecule_mol is None:
            raise ValueError(f"Invalid SMILES string: {molecule_smiles}")
        # return cast(bool, fragment_mol.HasSubstructMatch(molecule_mol))
        return cast(bool, molecule_mol.HasSubstructMatch(fragment_mol)) or cast(
            bool, molecule_mol.HasSubstructMatch(fragment_mol_withH)
        )

    def get_substructure_BBs(self, fragment: str) -> set[str]:
        """
        Get the set of building blocks (BBs) that the fragment matches.

        Args:
            fragment: SMILES string of the fragment.

        Returns:
            Set of building block SMILES strings that the fragment matches.
        """
        return set(bb for bb in self.BBs if self.is_strict_substructure(fragment, bb))
