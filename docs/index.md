# FragmentRetro: A Quadratic Retrosynthetic Method Based on Fragmentation Algorithms

FragmentRetro is a fragment-based retrosynthetic method with quadratic complexity, offering a scalable alternative to tree search and [DirectMultiStep](https://github.com/batistagroup/DirectMultiStep/tree/main).

## Quick Start

### Installation

You can install the package after cloning the repo:

```bash
uv venv --python 3.11.4
source .venv/bin/activate
uv pip install -e .
```

### Preprocess Stocks

To precompute properties and pattern fingerprints for faster screening, run the following commands:

```bash
python data/process_buyables.py
python data/process_stock.py
```

### Usage Example

Here's a quick example to obtain retrosynthesis solutions:

```python
from pathlib import Path
from fragmentretro.fragmenter import BRICSFragmenter  # or rBRICSFragmenter
from fragmentretro.retrosynthesis import Retrosynthesis
from fragmentretro.solutions import RetrosynthesisSolution

DATA_PATH = Path(__name__).parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
PRECOMPUTE_PATH = DATA_PATH / "precompute"
JSON_PRECOMPUTE_PATH = PRECOMPUTE_PATH / "n1_stock_properties.json"

target = "COc1ccc(F)c(-c2ccc(COc3cccc(C(CC(=O)O)C4CC4)c3)nc2CC(C)(C)C)c1"
fragmenter = BRICSFragmenter(target) # or rBRICSFragmenter
retro_tool = Retrosynthesis(fragmenter, original_BBs=None, mol_properties_path=JSON_PRECOMPUTE_PATH)
retro_tool.fragment_retrosynthesis()
retro_solution = RetrosynthesisSolution(retro_tool)
retro_solution.fill_solutions()
all_img = retro_solution.visualize_solutions(retro_solution.solutions)

# visualize first solution or simply run this in a notebook: all_img[0]
# with open("output.png", "wb") as f:
#     f.write(all_img[0].data)

# Or in a python script
# all_img[0].save("output.png", format="PNG")
```

### Graphical User Interface (GUI)

FragmentRetro's Graphical User Interface (GUI) is built using `ipywidget`. To access the GUI within a Jupyter notebook:

```python
from app.interface import display_gui

import logging
from fragmentretro.utils.logging_config import logger as fragment_logger
from app.logging_config import logger as app_logger

# Adjust the logging levels to control the verbosity of the logs or to suppress them
fragment_logger.setLevel(logging.WARNING)
app_logger.setLevel(logging.INFO) 

# To display the GUI without any initial SMILES input
app = display_gui()

# Alternatively, to display the GUI with a predefined SMILES string
app = display_gui(smiles="CCNCC")
```

## Development

We welcome any contributions, feel free to clone the repo and create a PR. We recommend using [uv](https://docs.astral.sh/uv/getting-started/installation/):

```bash
uv venv --python 3.11.4
source .venv/bin/activate
uv pip install -e ".[dev]"
```

## License

- Code: MIT License
- Paper content ([arXiv preprint](https://arxiv.org/abs/xxxxx)): CC-BY 4.0
