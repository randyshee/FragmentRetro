# Helper Functions

This module provides helper functions for working with SMILES strings.

## Example Use

`canonicalize_smiles` is for canonicalizing SMILES strings.

```python
from fragmentretro.utils.helpers import canonicalize_smiles

smiles = 'C[C@@H](O)C(=O)O'
canonical_smiles = canonicalize_smiles(smiles)
print(f"Original SMILES: {smiles}")
print(f"Canonical SMILES: {canonical_smiles}")

```

`replace_dummy_atoms_regex` is for replacing dummy atoms with hydrogen atoms for pattern fingerprint screening. See how it's used in the `filter_compounds` function in the `CompoundFilter` class.

```python
from fragmentretro.utils.helpers import replace_dummy_atoms_regex

smiles_with_dummy = '[5*]N1CCC[C@@]1([13*])C'
smiles_without_dummy = replace_dummy_atoms_regex(smiles_with_dummy)
print(f"SMILES with dummy atoms: {smiles_with_dummy}")
print(f"SMILES without dummy atoms: {smiles_without_dummy}")
```

`remove_indices_before_dummy` is for removing indices before dummy atoms. This is to record processed fragment SMILES strings in the most general format. See how it's used in the `Retrosythesis` class as well.

```python
from fragmentretro.utils.helpers import remove_indices_before_dummy

smiles_with_indices = '[5*]N1CCC[C@@]1([13*])C'
smiles_without_indices = remove_indices_before_dummy(smiles_with_indices)
print(f"SMILES with indices: {smiles_with_indices}")
print(f"SMILES without indices: {smiles_without_indices}")
```

## Source Code

::: FragmentRetro.utils.helpers
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - canonicalize_smiles
        - count_heavy_atoms
        - sort_by_heavy_atoms
        - replace_dummy_atoms_regex
        - remove_indices_before_dummy
