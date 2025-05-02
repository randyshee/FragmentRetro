from IPython.display import display

from app.gui.controller import GuiController
from app.gui.layout import gui_layout
from app.gui.state import AppState
from app.gui.widgets import (
    fragment_comb_dropdown,
    image_display_area,
    next_smiles_button,
    output_area,
    prev_smiles_button,
    smiles_display_area,
    smiles_pagination_label,
    solution_dropdown,
    solution_output_area,
    sort_smiles_button,
    target_smiles_input,
)


# --- Main Function to Display GUI ---
def display_gui(smiles: str | None = None) -> None:
    """Instantiates state, controller, connects handlers, and displays the GUI."""

    # 1. Instantiate State
    app_state = AppState()

    # 2. Instantiate Controller (pass state)
    # The controller __init__ already references the imported widget instances
    controller = GuiController(app_state)

    # 3. Register Event Handlers
    controller.register_event_handlers()

    # 4. Initial UI Setup / Reset
    if smiles:
        target_smiles_input.value = smiles

    # Ensure initial state is clean (using widget instances directly)
    solution_dropdown.options = []
    solution_dropdown.value = None
    solution_dropdown.disabled = True
    image_display_area.clear_output(wait=False)
    solution_output_area.clear_output(wait=False)
    output_area.clear_output(wait=False)

    fragment_comb_dropdown.options = []
    fragment_comb_dropdown.value = None
    fragment_comb_dropdown.disabled = True
    prev_smiles_button.disabled = True
    next_smiles_button.disabled = True
    smiles_pagination_label.value = "0 of 0"
    sort_smiles_button.disabled = True
    smiles_display_area.clear_output(wait=False)

    # Reset backend state as well
    app_state.reset_run_state()

    # 5. Display Layout
    display(gui_layout)  # type: ignore
