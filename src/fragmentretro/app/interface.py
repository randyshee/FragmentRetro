from app.gui.controller import GuiController
from app.gui.layout import gui_layout
from app.gui.state import AppState
from app.gui.widgets import target_smiles_input
from IPython.display import display


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

    controller.reset_ui_outputs()  # Reset UI elements via controller
    app_state.reset_run_state()  # Reset backend state as well

    # 5. Display Layout
    display(gui_layout)  # type: ignore
