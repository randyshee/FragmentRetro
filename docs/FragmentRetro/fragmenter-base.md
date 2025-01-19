# Base Class for Fragmenter

## Example Use

You will only have to define the abstract methods to find fragmentation bonds and to break these bonds. Other methods like building fragment graph (`networkx` graph), visualization, and how to get SMILES string given a combination of fragments are provided in this base class.

```python
from FragmentRetro.fragmenter_base import Fragmenter

# suppose you have a new way to find and break bonds: `FindTestBonds` and `BreakTestBonds`

class TestFragmenter(Fragmenter):
    def __init__(self, smiles: str) -> None:
        """
        Initialize with SMILES string.

        Args:
            smiles: SMILES string of molecule to fragment
        """
        super().__init__(smiles)

    def _find_fragmentation_bonds(self, mol: Mol) -> BondsType:
        return FindTestBonds(mol)

    def _break_bonds(self, mol: Mol, bonds: BondsType) -> Mol:
        return BreakTestBonds(mol, bonds)
```

## Source Code

::: FragmentRetro.fragmenter_base
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - Fragmenter