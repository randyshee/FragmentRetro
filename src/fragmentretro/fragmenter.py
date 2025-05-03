"""Module for fragmenting molecules using BRICS (add other algorithms later)."""

from fragmentation.rBRICS import BreakrBRICSBonds, FindrBRICSBonds
from rdkit.Chem import Mol
from rdkit.Chem.BRICS import BreakBRICSBonds, FindBRICSBonds

from FragmentRetro.fragmenter_base import Fragmenter
from FragmentRetro.utils.type_definitions import BondType


class BRICSFragmenter(Fragmenter):
    def __init__(self, smiles: str) -> None:
        """
        Initialize with SMILES string.

        Args:
            smiles: SMILES string of molecule to fragment
        """
        super().__init__(smiles)

    def _find_fragmentation_bonds(self, mol: Mol) -> list[BondType]:
        return list(FindBRICSBonds(mol))

    def _break_bonds(self, mol: Mol, bonds: list[BondType]) -> Mol:
        return BreakBRICSBonds(mol, bonds)


class rBRICSFragmenter(Fragmenter):
    def __init__(self, smiles: str) -> None:
        """
        Initialize with SMILES string.

        Args:
            smiles: SMILES string of molecule to fragment
        """
        super().__init__(smiles)

    def _find_fragmentation_bonds(self, mol: Mol) -> list[BondType]:
        return list(FindrBRICSBonds(mol))  # type: ignore[no-untyped-call]

    def _break_bonds(self, mol: Mol, bonds: list[BondType]) -> Mol:
        return BreakrBRICSBonds(mol, bonds)  # type: ignore[no-untyped-call]
