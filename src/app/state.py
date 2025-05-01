from PIL.Image import Image as PILImage  # For type hinting

from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution
from FragmentRetro.utils.type_definitions import CombType, SolutionType


class AppState:
    """Holds the state of the FragmentRetro GUI application."""

    def __init__(self) -> None:
        self.retro_solution: RetrosynthesisSolution | None = None
        self.valid_images_cache: list[PILImage] = []
        self.retro_tool: Retrosynthesis | None = None
        self.selected_fragment_comb: CombType | None = None  # Assuming fragment comb is tuple of strings
        self.current_smiles_list: list[str] = []
        self.current_smiles_index: int = 0
        self.displayable_solutions: list[
            SolutionType
        ] = []  # Type hint might need refinement based on solution structure
        self.is_smiles_sorted: bool = False  # NEW: Track if the current SMILES list is sorted

    def reset_run_state(self) -> None:
        """Resets state related to a specific retrosynthesis run."""
        self.retro_solution = None
        self.valid_images_cache = []
        self.retro_tool = None
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0
        self.displayable_solutions = []
        self.is_smiles_sorted = False  # Reset sorting state

    def reset_display_state(self) -> None:
        """Resets state related to the solution display area."""
        self.valid_images_cache = []
        self.displayable_solutions = []
        # Keep retro_solution and retro_tool, but reset downstream state
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0
        self.is_smiles_sorted = False  # Reset sorting state

    def reset_smiles_viewer_state(self) -> None:
        """Resets state specific to the SMILES viewer section."""
        self.selected_fragment_comb = None
        self.current_smiles_list = []
        self.current_smiles_index = 0
        self.is_smiles_sorted = False  # Reset sorting state
