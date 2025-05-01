from typing import Any, Optional, Tuple

from PIL.Image import Image as PILImage  # For type hinting

from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution


class AppState:
    """Holds the state of the FragmentRetro GUI application."""

    def __init__(self):
        self.retro_solution: Optional[RetrosynthesisSolution] = None
        self.valid_images_cache: list[PILImage] = []
        self.retro_tool: Optional[Retrosynthesis] = None
        self.selected_fragment_comb: Optional[Tuple[str, ...]] = None  # Assuming fragment comb is tuple of strings
        self.current_smiles_list: list[str] = []
        self.current_smiles_index: int = 0
        self.displayable_solutions: list[Any] = []  # Type hint might need refinement based on solution structure

    def reset_run_state(self):
        """Resets state related to a specific retrosynthesis run."""
        self.retro_solution = None
        self.valid_images_cache = []
        self.retro_tool = None
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0
        self.displayable_solutions = []

    def reset_display_state(self):
        """Resets state related to the solution display area."""
        self.valid_images_cache = []
        self.displayable_solutions = []
        # Keep retro_solution and retro_tool, but reset downstream state
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0

    def reset_smiles_viewer_state(self):
        """Resets state specific to the SMILES viewer section."""
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0
