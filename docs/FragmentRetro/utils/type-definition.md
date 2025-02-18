# Type Definitions

[type_definitions.py](/src/FragmentRetro/type_definitions.py) is a good place to store your type definitions to make code clean and pass mypy tests.

## Example Usage

```python
def _find_fragmentation_bonds(mol: Mol) -> list[BondType]:
    return list(FindBRICSBonds(mol))
```

where `BondType` is defined in [type_definitions.py](/src/FragmentRetro/utils/type_definitions.py) as:

```python
BondType: TypeAlias = tuple[tuple[int, int], tuple[str, str]]
```

## Source Code

::: FragmentRetro.utils.type_definitions
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
