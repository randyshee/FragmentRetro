# Substructure Match

Check if a fragment (which can contain dummy or any atoms represented by "*") is a strict substructure of a molecule. A strict substructure means that no additional atoms should be present in the molecule other than those explicitly defined in the fragment. Extra atoms (including hydrogen atoms) are only allowed at the sites of dummy atoms in the fragment.

## Example Use

```python
from FragmentRetro.substructure_matcher.SubstructureMatcher import is_strict_substructure

fragment_smiles = "[4*]CCN[5*]",
molecule_smiles = "CCCNC"

print(is_strict_substructure(fragment_smiles, molecule_smiles))
# should print True
```

## Source Code

::: FragmentRetro.substructure_matcher
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - SubstructureMatcher
