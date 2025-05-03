# Type Definitions

[typing.py](/src/fragmentretro/typing.py) is a good place to store your type definitions to make code clean and pass mypy tests.

## Example Usage

```python
def _find_fragmentation_bonds(mol: Mol) -> list[BondType]:
    return list(FindBRICSBonds(mol))
```

where `BondType` is defined in [typing.py](/src/fragmentretro/typing.py) as:

```python
BondType: TypeAlias = tuple[tuple[int, int], tuple[str, str]]
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
