# Derived Class for Fragmenter

Once you have the fragmentation algorithms such as BRICS, you can get derived class from Fragmenter base class.

## Example Use

```python
from fragmentretro.fragmenter import BRICSFragmenter

smiles = "COc1ccc(-n2nccn2)c(C(=O)N2CCC[C@@]2(C)c2nc3c(C)c(Cl)ccc3[nH]2)c1"
fragmenter = BRICSFragmenter(smiles)

# Custom figure size
fragmenter.visualize(figsize=(20, 10), verbose=False, with_indices=True)
```

## Source Code

::: FragmentRetro.fragmenter
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - BRICSFragmenter
        - rBRICSFragmenter
