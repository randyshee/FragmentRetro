"""Module for fragmenting molecules using BRICS (add other algorithms later)."""

from rdkit.Chem import Mol
from rdkit.Chem.BRICS import BreakBRICSBonds, FindBRICSBonds

from fragmentretro.fragmentation.r_brics import break_r_brics_bonds, find_brics_bonds
from fragmentretro.fragmenter_base import Fragmenter
from fragmentretro.typing import BondType


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
        return list(find_brics_bonds(mol))  # type: ignore[no-untyped-call]

    def _break_bonds(self, mol: Mol, bonds: list[BondType]) -> Mol:
        return break_r_brics_bonds(mol, bonds)  # type: ignore[no-untyped-call]
