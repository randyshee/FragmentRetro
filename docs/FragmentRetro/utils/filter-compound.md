# Compound Filtering

This module provides functions for filtering compounds based on precomputed molecular properties (number of heavy atoms, number of rings, and pattern fingerprints).

## Example Use

The `precompute_properties` function computes properties for a large set of compounds and saves them to a JSON file.  The example below shows how to use it with a subset of the `n1-stock.txt` file.

```python
from pathlib import Path
from FragmentRetro.utils.filter_compound import precompute_properties

DATA_PATH = Path(__file__).parent.parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
PRECOMPUTE_PATH = DATA_PATH / "precompute"
MOL_PROPERTIES_PATH = PRECOMPUTE_PATH / "n1_stock_properties_subset.json"

# Create the directory if it doesn't exist
PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)

with open(PAROUTES_PATH / "n1-stock.txt") as f:
    n1_stock = [line.strip() for line in f]
n1_stock_subset = n1_stock[:500]
precompute_properties(n1_stock_subset, MOL_PROPERTIES_PATH, fpSize=2048)
```

The `CompoundFilter` class is used to initialize a filter, and the `get_filtered_BBs` method retrieves the screened building blocks. The following example demonstrates filtering using a specific fragment SMILES string.

```python
from FragmentRetro.utils.filter_compound import CompoundFilter
from FragmentRetro.utils.helpers import replace_dummy_atoms_regex

fragment_smiles_list = ["[5*]N1CCC[C@@]1([13*])C", "[4*]CCN[5*]", "[4*]C[8*]", "[*]C[*]", "[3*]O[3*]"]
fragment_smiles = fragment_smiles_list[1]

compound_filter = CompoundFilter(MOL_PROPERTIES_PATH, fpSize=2048)
filtered_indices, filtered_BBs = compound_filter.get_filtered_BBs(fragment_smiles)
```

## Source Code

::: FragmentRetro.utils.filter_compound
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - get_mol_properties
        - precompute_properties
        - CompoundFilter
