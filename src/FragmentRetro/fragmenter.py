"""Module for fragmenting molecules using BRICS (add other algorithms later)."""

from rdkit.Chem import Mol
from rdkit.Chem.BRICS import BreakBRICSBonds, FindBRICSBonds

from .fragmenter_base import Fragmenter
from .type_definitions import BondsType


class BRICSFragmenter(Fragmenter):
    def __init__(self, smiles: str) -> None:
        """
        Initialize with SMILES string.

        Args:
            smiles: SMILES string of molecule to fragment
        """
        super().__init__(smiles)

    def _find_fragmentation_bonds(self, mol: Mol) -> BondsType:
        return list(FindBRICSBonds(mol))

    def _break_bonds(self, mol: Mol, bonds: BondsType) -> Mol:
        return BreakBRICSBonds(mol, bonds)
