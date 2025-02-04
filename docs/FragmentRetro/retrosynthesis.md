# Retrosynthesis

The `Retrosynthesis` class identifies all possible retrosynthesis solutions for a given molecule (after the valid fragment combinations are passed to the `RetrosynthesisSolution` class). It uses a fragmenter to break down the molecule into fragments, then recombines these fragments in every possible way, checking each combination against a set of allowed building blocks.

## Example Use

```python
from FragmentRetro.fragmenter import BRICSFragmenter
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution

# Define your building block set (SMILES strings)
original_BBs = set(['Brc1cc(OC)ccc1-n1nccn1',
                    'BrC(Br)=O',
                    'BrN1CCC[C@@]1(Br)C',
                    'Brc1nc2c(C)c(Cl)ccc2[nH]1'])

# Initialize a fragmenter with the target SMILES
target = "COc1ccc(-n2nccn2)c(C(=O)N2CCC[C@@]2(C)c2nc3c(C)c(Cl)ccc3[nH]2)c1"
fragmenter = BRICSFragmenter(target)

# Initialize the retrosynthesis tool with the fragmenter and building blocks
retro_tool = Retrosynthesis(fragmenter, original_BBs)

# Run the retrosynthesis algorithm
retro_tool.fragment_retrosynthesis()

# Pass the retrosynthesis tool to the solution tool to process the results
retro_solution = RetrosynthesisSolution(retro_tool)
retro_solution.fill_solutions()
print(retro_solution.solutions)
# Expected output:
# [[(0, 1, 2), (3,), (4,), (5,)],
#  [(0, 1), (2,), (3,), (4,), (5,)],
#  [(0,), (1, 2), (3,), (4,), (5,)],
#  [(0,), (1,), (2,), (3,), (4,), (5,)]]

# Visualize the first solution
retro_solution.visualize_solutions(retro_solution.solutions, molsPerRow=4)[0]

# Get the valid building blocks for a fragment combination 
# retro_tool.comb_bbs_dict[(0, 1 ,2)]
```

## Source Code

::: FragmentRetro.retrosynthesis
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - Retrosynthesis
