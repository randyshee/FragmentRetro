# Retrosynthesis Solution

The `RetrosynthesisSolution` class is designed to find and visualize all possible retrosynthesis solutions given a set of fragmented molecules and their valid combinations. `Retrosynthesis` class is used to get valid combinations of fragments and `Fragmenter` class is used to obtain the fragment SMILES for visualization.

## Example Use

```python
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.fragmenter import BRICSFragmenter
from FragmentRetro.solutions import RetrosynthesisSolution

# Example SMILES
smiles = "COc1ccc(-n2nccn2)c(C(=O)N2CCC[C@@]2(C)c2nc3c(C)c(Cl)ccc3[nH]2)c1"
fragmenter = BRICSFragmenter(smiles)

retro_tool = Retrosynthesis(fragmenter=fragmenter, original_BBs=set())

# no need to run retrosynthesis for this example usage
# will give an example solution in the next few lines
# retro_tool.fragment_retrosynthesis 

# Example Solution
retro_solution = RetrosynthesisSolution(retro_tool)
retro_solution.solutions = [[[0, 1, 2], [3], [4], [5]]]

# Visualize solutions
images = retro_solution.visualize_solutions(retro_solution.solutions, molsPerRow=4)
for i, img in enumerate(images):
    img.save(f"solution_{i}.png")
    print(f"Solution {i} saved to solution_{i}.png")
```

## Source Code

::: FragmentRetro.solutions
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - RetrosynthesisSolution
