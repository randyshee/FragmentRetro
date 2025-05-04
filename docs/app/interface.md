# FragmentRetro Interface

This is the GUI Implementation for FragmentRetro.

## Example Use

In a Jupyter notebook:

```python
from fragmentretro.app.interface import display_gui

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
