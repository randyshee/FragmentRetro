# Type Definitions

`fragmentretro.typing` is a good place to store your type definitions to make code clean and pass mypy tests.

## Example Usage

```python
from fragmentretro.typing import BondType

def _find_fragmentation_bonds(mol: Mol) -> list[BondType]:
    return list(FindBRICSBonds(mol))
```

where `BondType` is defined in `fragmentretro.typing` as:

```python
type BondType = tuple[tuple[int, int], tuple[str, str]]
```

## Source Code

::: fragmentretro.typing
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - BondType
        - AtomMappingType
        - CombType
        - SolutionType
        - BBsType
        - StageCombDictType
        - CombBBsDictType
        - FragmentBBsDictType
        - FilterIndicesType
        - CombFilterIndicesDictType
        - MolProperties
