import contextlib
import io
from pathlib import Path
from typing import cast

import ipywidgets as widgets
from IPython.display import display
from PIL.Image import Image as PILImage
from rdkit import Chem
from rdkit.Chem import Draw

from fragmentretro.app.gui.state import AppState
from fragmentretro.app.gui.widgets import (
    core_factor_input,
    display_button,
    file_path_input,
    filter_checkbox,
    fragment_comb_dropdown,
    fragment_count_input,
    fragmenter_choice,
    image_display_area,
    next_smiles_button,
    num_cores_input,
    output_area,
    parallelize_checkbox,
    prev_smiles_button,
    run_button,
    smiles_display_area,
    smiles_pagination_label,
    solution_dropdown,
    solution_output_area,
    sort_smiles_button,
    target_smiles_input,
)
from fragmentretro.fragmenter import BRICSFragmenter, rBRICSFragmenter
from fragmentretro.fragmenter_base import Fragmenter
from fragmentretro.retrosynthesis import Retrosynthesis
from fragmentretro.solutions import RetrosynthesisSolution
from fragmentretro.typing import CombType, SolutionType
from fragmentretro.utils.helpers import sort_by_heavy_atoms
from fragmentretro.utils.logging_config import logger


class GuiController:
    """Manages GUI event handling and interactions."""

    def __init__(self, state: AppState):
        self.state = state
        self.target_smiles_input = target_smiles_input
        self.fragmenter_choice = fragmenter_choice
        self.file_path_input = file_path_input
        self.parallelize_checkbox = parallelize_checkbox
        self.num_cores_input = num_cores_input
        self.core_factor_input = core_factor_input
        self.run_button = run_button
        self.output_area = output_area
        self.filter_checkbox = filter_checkbox
        self.fragment_count_input = fragment_count_input
        self.display_button = display_button
        self.solution_output_area = solution_output_area
        self.solution_dropdown = solution_dropdown
        self.image_display_area = image_display_area
        self.fragment_comb_dropdown = fragment_comb_dropdown
        self.smiles_display_area = smiles_display_area
        self.prev_smiles_button = prev_smiles_button
        self.smiles_pagination_label = smiles_pagination_label
        self.next_smiles_button = next_smiles_button
        self.sort_smiles_button = sort_smiles_button

    def register_event_handlers(self) -> None:
        """Connects widget events to controller methods."""
        self.parallelize_checkbox.observe(self.handle_parallelize_change, names="value")
        self.filter_checkbox.observe(self.handle_filter_change, names="value")
        self.run_button.on_click(self.run_retrosynthesis_on_click)
        self.fragment_comb_dropdown.observe(self.on_fragment_comb_select, names="value")
        self.prev_smiles_button.on_click(self.on_prev_smiles_click)
        self.next_smiles_button.on_click(self.on_next_smiles_click)
        self.sort_smiles_button.on_click(self.on_sort_smiles_click)
        self.display_button.on_click(self.display_solutions_on_click)
        self.solution_dropdown.observe(self.on_solution_select, names="value")

    def handle_parallelize_change(self, change: dict[str, bool]) -> None:
        """Callback to enable/disable core inputs based on parallelize checkbox."""
        is_parallel = change.get("new", False)
        self.num_cores_input.disabled = not is_parallel
        self.core_factor_input.disabled = not is_parallel

    def handle_filter_change(self, change: dict[str, bool]) -> None:
        """Callback to enable/disable fragment count input based on filter checkbox."""
        is_filtered = change.get("new", False)
        self.fragment_count_input.disabled = not is_filtered
        if not is_filtered:
            self.fragment_count_input.value = None  # type: ignore

    def reset_ui_outputs(self) -> None:
        """Resets output areas, dropdowns, and buttons to their initial states."""
        self.solution_dropdown.options = []
        self.solution_dropdown.value = None
        self.solution_dropdown.disabled = True
        self.image_display_area.clear_output(wait=False)
        self.solution_output_area.clear_output(wait=False)

        self.fragment_comb_dropdown.options = []
        self.fragment_comb_dropdown.value = None
        self.fragment_comb_dropdown.disabled = True
        self.prev_smiles_button.disabled = True
        self.next_smiles_button.disabled = True
        self.smiles_pagination_label.value = "0 of 0"
        self.sort_smiles_button.disabled = True
        self.smiles_display_area.clear_output(wait=False)
        with self.smiles_display_area:
            logger.info("[GUI] Perform an action (Run/Display) to populate this area.")

    def run_retrosynthesis_on_click(self, b: widgets.Button) -> None:
        """Handles the click event for the 'Run Retrosynthesis' button."""
        self.state.reset_run_state()
        self.reset_ui_outputs()  # Reset UI elements

        self.output_area.clear_output(wait=True)  # Clear main output specifically for run
        with self.output_area:
            logger.info("[GUI] Starting retrosynthesis...")

        target: str = self.target_smiles_input.value
        fragmenter_name: str = self.fragmenter_choice.value
        json_path_str: str = self.file_path_input.value
        parallelize: bool = self.parallelize_checkbox.value
        num_cores: int | None = self.num_cores_input.value  # type: ignore
        core_factor: int = self.core_factor_input.value  # type: ignore

        if not target:
            with self.output_area:
                logger.error("[GUI] ERROR: Target SMILES cannot be empty.")
            return

        if json_path_str:
            mol_properties_path = Path(json_path_str)
            if not mol_properties_path.is_file():
                with self.output_area:
                    logger.error(f"[GUI] ERROR: Properties file not found at {mol_properties_path}")
                return
        try:
            logger.info("[GUI] Running fragmentation...")
            fragmenter: Fragmenter | None = None
            if fragmenter_name == "BRICSFragmenter":
                fragmenter = BRICSFragmenter(target)
            elif fragmenter_name == "rBRICSFragmenter":
                fragmenter = rBRICSFragmenter(target)
            else:
                with self.output_area:
                    logger.error(f"[GUI] ERROR: Unknown fragmenter type '{fragmenter_name}'")
                return

            with self.output_area:
                logger.info(f"[GUI] Using Fragmenter: {fragmenter_name}")
                logger.info(f"[GUI] Using Properties: {mol_properties_path}")
                logger.info(f"[GUI] Parallelize: {parallelize}, Num Cores: {num_cores}, Core Factor: {core_factor}")

            retro_tool = Retrosynthesis(
                fragmenter,
                mol_properties_path=mol_properties_path,
                parallelize=parallelize,
                num_cores=num_cores,
                core_factor=core_factor,
            )
            with self.output_area:
                logger.info("[GUI] Running retrosynthesis...")
                retro_tool.fragment_retrosynthesis()
                retro_solution = RetrosynthesisSolution(retro_tool)
                retro_solution.fill_solutions()
                logger.info(f"[GUI] Found {len(retro_solution.solutions)} solution(s).")

            self.state.retro_solution = retro_solution
            self.state.retro_tool = retro_tool
            with self.output_area:
                logger.info("[GUI] Retrosynthesis complete. Ready to display solutions and browse fragment SMILES.")

        except Exception as e:
            with self.output_area:
                logger.error(f"[GUI] ERROR during retrosynthesis: {e}", exc_info=True)

    def update_fragment_comb_dropdown(self, solution: SolutionType | None) -> None:
        """Populates the fragment comb dropdown based on a single solution."""
        with contextlib.suppress(ValueError):
            self.fragment_comb_dropdown.unobserve(self.on_fragment_comb_select, names="value")

        if solution:
            self.fragment_comb_dropdown.options = [(str(comb), comb) for comb in solution]
            self.fragment_comb_dropdown.disabled = False
            self.smiles_display_area.clear_output(wait=False)
            logger.debug("[GUI] Select a fragment combination to view SMILES.")
        else:
            self.fragment_comb_dropdown.options = []
            self.fragment_comb_dropdown.disabled = True
            self.smiles_display_area.clear_output(wait=False)
            logger.debug("[GUI] No fragment combinations in the selected solution.")

        self.fragment_comb_dropdown.value = None
        self.prev_smiles_button.disabled = True
        self.next_smiles_button.disabled = True
        self.smiles_pagination_label.value = "0 of 0"
        self.state.selected_fragment_comb = None
        self.state.current_smiles_list = []
        self.state.current_smiles_index = 0

        self.fragment_comb_dropdown.observe(self.on_fragment_comb_select, names="value")

    def update_smiles_display(self) -> None:
        """Updates the SMILES display area and pagination controls."""
        smiles_list: list[str] = self.state.current_smiles_list
        index: int = self.state.current_smiles_index
        total_smiles: int = len(smiles_list)

        with self.smiles_display_area:
            self.smiles_display_area.clear_output(wait=True)
            if not smiles_list:
                logger.info("[GUI] No SMILES found for this combination.")
                self.smiles_pagination_label.value = "0 of 0"
                self.prev_smiles_button.disabled = True
                self.next_smiles_button.disabled = True
                return

            current_smiles: str = smiles_list[index]
            try:
                mol = Chem.MolFromSmiles(current_smiles)
                if mol:
                    img: PILImage = Draw.MolToImage(mol, size=(300, 300))
                    bio = io.BytesIO()
                    img.save(bio, format="PNG")
                    image_widget = widgets.Image(value=bio.getvalue(), format="png", width=300, height=300)
                    display(image_widget)  # type: ignore
                else:
                    logger.warning(f"[GUI] Invalid SMILES: {current_smiles}")
            except Exception as e:
                logger.error(f"[GUI] Error generating image for SMILES {current_smiles}: {e}")

        self.smiles_pagination_label.value = f"{index + 1} of {total_smiles}"
        self.prev_smiles_button.disabled = index == 0
        self.next_smiles_button.disabled = index >= total_smiles - 1

    def on_fragment_comb_select(self, change: dict[str, str | CombType]) -> None:
        """Handles selection changes in the fragment combination dropdown."""
        if change.get("type") == "change" and change.get("name") == "value":
            selected_comb: CombType = cast(CombType, change.get("new"))
            retro_tool: Retrosynthesis | None = self.state.retro_tool

            if selected_comb is None or retro_tool is None or not hasattr(retro_tool, "comb_bbs_dict"):
                self.state.reset_smiles_viewer_state()
                self.sort_smiles_button.disabled = True
            else:
                smiles_data: set[str] = retro_tool.comb_bbs_dict.get(selected_comb, set())
                smiles_list: list[str] = list(smiles_data)
                self.state.selected_fragment_comb = selected_comb
                self.state.current_smiles_list = smiles_list
                self.state.current_smiles_index = 0
                self.state.is_smiles_sorted = False
                self.sort_smiles_button.disabled = not smiles_list

            self.update_smiles_display()

    def on_prev_smiles_click(self, b: widgets.Button) -> None:
        """Handles clicks on the 'Previous' SMILES button."""
        if self.state.current_smiles_index > 0:
            self.state.current_smiles_index -= 1
            self.update_smiles_display()

    def on_next_smiles_click(self, b: widgets.Button) -> None:
        """Handles clicks on the 'Next' SMILES button."""
        smiles_list: list[str] = self.state.current_smiles_list
        if self.state.current_smiles_index < len(smiles_list) - 1:
            self.state.current_smiles_index += 1
            self.update_smiles_display()

    def on_sort_smiles_click(self, b: widgets.Button) -> None:
        """Handles clicks on the 'Sort by Heavy Atoms' button."""
        if not self.state.is_smiles_sorted and self.state.current_smiles_list:
            self.state.current_smiles_list = sort_by_heavy_atoms(self.state.current_smiles_list)
            self.state.is_smiles_sorted = True
            self.state.current_smiles_index = 0
            self.update_smiles_display()
            self.sort_smiles_button.disabled = True

    def display_solutions_on_click(self, b: widgets.Button) -> None:
        """Handles the click event for the 'Display Solutions' button."""
        self.reset_ui_outputs()  # Reset UI elements
        self.state.reset_display_state()  # Reset backend display state
        with self.solution_output_area:  # Log specifically to solution output area after clear
            logger.debug("[GUI] Attempting to display solutions...")

        retro_solution: RetrosynthesisSolution | None = self.state.retro_solution
        use_filter: bool = self.filter_checkbox.value
        fragment_count: int | None = self.fragment_count_input.value if use_filter else None  # type: ignore

        if retro_solution is None:
            msg = "No retrosynthesis results available. Please run retrosynthesis first."
            with self.solution_output_area:
                logger.warning(f"[GUI] {msg}")
            return

        filtered_solutions: list[SolutionType] = []
        if use_filter and (fragment_count is None or not isinstance(fragment_count, int) or fragment_count <= 0):
            msg = "Filter checkbox is checked, but fragment count is invalid. Displaying all solutions."
            with self.solution_output_area:
                logger.warning(f"[GUI] {msg}")
            filtered_solutions = retro_solution.solutions
        elif use_filter:
            filtered_solutions = [sol for sol in retro_solution.solutions if len(sol) == fragment_count]
            msg = f"Filtering for solutions with exactly {fragment_count} fragments."
            with self.solution_output_area:
                logger.info(f"[GUI] {msg}")
        else:
            msg = "Filter checkbox is unchecked. Displaying all available solutions."
            with self.solution_output_area:
                logger.info(f"[GUI] {msg}")
            filtered_solutions = retro_solution.solutions

        if not filtered_solutions:
            msg = f"No solutions found matching the criteria (count: {fragment_count})."
            with self.solution_output_area:
                logger.info(f"[GUI] {msg}")
            return

        with self.solution_output_area:
            logger.info(f"[GUI] Visualizing {len(filtered_solutions)} solution(s)...")

        try:
            solution_images: list[PILImage] = retro_solution.visualize_solutions(filtered_solutions)
            if not solution_images:
                msg = "Visualization did not produce any images."
                with self.solution_output_area:
                    logger.info(f"[GUI] {msg}")
                return

            valid_images: list[PILImage] = []
            displayable_solutions: list[SolutionType] = []
            original_indices: list[int] = []
            for i, img in enumerate(solution_images):
                if img:
                    valid_images.append(img)
                    current_filtered_solution = filtered_solutions[i]
                    displayable_solutions.append(current_filtered_solution)
                    try:
                        original_solution_index = retro_solution.solutions.index(current_filtered_solution)
                        original_indices.append(original_solution_index)
                    except ValueError:
                        with self.solution_output_area:
                            logger.error(f"[GUI] Could not find filtered solution {i} in original solutions list.")
                        original_indices.append(-1)
                else:
                    with self.solution_output_area:
                        logger.warning(
                            f"[GUI] Null image returned by visualize_solutions for filtered solution index {i}"
                        )

            num_valid_images = len(valid_images)
            self.state.valid_images_cache = valid_images
            self.state.displayable_solutions = displayable_solutions

            if num_valid_images == 0:
                msg = "Visualization did not produce any valid images."
                with self.solution_output_area:
                    logger.info(f"[GUI] {msg}")
                return

            with self.solution_output_area:
                logger.info(f"[GUI] Generated {num_valid_images} image(s). Use dropdown to view.")

            dropdown_options: list[tuple[str, int]] = [
                (f"Solution {original_indices[i] + 1}", i) for i in range(num_valid_images) if original_indices[i] != -1
            ]
            self.solution_dropdown.options = dropdown_options
            self.solution_dropdown.value = 0 if dropdown_options else None
            self.solution_dropdown.disabled = not dropdown_options
            self.solution_dropdown.observe(self.on_solution_select, names="value")

            if num_valid_images > 0 and self.state.displayable_solutions:
                initial_solution: SolutionType = self.state.displayable_solutions[0]
                self.update_fragment_comb_dropdown(initial_solution)
            else:
                self.update_fragment_comb_dropdown(None)

            if valid_images:
                self.image_display_area.clear_output(wait=True)
                with self.image_display_area:
                    display(valid_images[0])  # type: ignore

            with self.solution_output_area:
                logger.info("[GUI] Solutions displayed.")

        except Exception as e:
            with self.solution_output_area:
                logger.error(f"[GUI] An error occurred during solution visualization: {e}", exc_info=True)

    def on_solution_select(self, change: dict[str, str | int]) -> None:
        """Callback function for solution dropdown selection changes."""
        if change.get("type") == "change" and change.get("name") == "value":
            selected_index: int = cast(int, change.get("new"))
            valid_images_cache: list[PILImage] = self.state.valid_images_cache
            displayable_solutions: list[SolutionType] = self.state.displayable_solutions

            if selected_index is not None and 0 <= selected_index < len(valid_images_cache):
                self.image_display_area.clear_output(wait=True)
                with self.image_display_area:
                    display(valid_images_cache[selected_index])  # type: ignore

                if selected_index < len(displayable_solutions):
                    selected_solution: SolutionType = displayable_solutions[selected_index]
                    self.update_fragment_comb_dropdown(selected_solution)
                    self.state.reset_smiles_viewer_state()
                else:
                    with self.solution_output_area:
                        logger.warning(
                            f"[GUI] Selected index {selected_index} out of bounds for displayable_solutions."
                        )
                    self.update_fragment_comb_dropdown(None)
                    self.state.reset_smiles_viewer_state()

            elif selected_index is None:
                self.image_display_area.clear_output(wait=True)
                with self.solution_output_area:
                    logger.info("[GUI] No solution selected or available.")
                self.update_fragment_comb_dropdown(None)
                self.state.reset_smiles_viewer_state()
                self.sort_smiles_button.disabled = True
