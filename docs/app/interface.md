# FragmentRetro Interface

This is the GUI Implementation for FragmentRetro.

## Example Use

In a Jupyter notebook:

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

## Source Code

::: fragmentretro.app.interface
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - display_gui

::: fragmentretro.app.gui.controller
    handler: python
    options:
      show_root_heading: true
      show_source: true
      members:
        - GuiController
