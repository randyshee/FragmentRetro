import io  # Add io for image handling
from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output, display
from rdkit import Chem  # Add RDKit Chem
from rdkit.Chem import Draw  # Add RDKit Draw

from app.logging_config import logger
from FragmentRetro.fragmenter import BRICSFragmenter, rBRICSFragmenter
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution
from FragmentRetro.utils.helpers import sort_by_heavy_atoms

# --- Widgets ---

# Input fields
target_smiles_input = widgets.Text(
    placeholder="Enter target molecule SMILES string",
    description="Target:",
    layout=widgets.Layout(width="90%"),
    style={"description_width": "initial"},
)

fragmenter_choice = widgets.Dropdown(
    options=["BRICSFragmenter", "rBRICSFragmenter"],
    value="BRICSFragmenter",
    description="Fragmenter:",
    style={"description_width": "initial"},
)

file_path_input = widgets.Text(
    value="data/precompute/n1_stock_properties.json",  # Example default path relative to notebook location potentially
    placeholder="Path to BBs or properties JSON",
    description="Properties JSON:",
    layout=widgets.Layout(width="90%"),
    style={"description_width": "initial"},
)

# Parallelization options
parallelize_checkbox = widgets.Checkbox(value=False, description="Parallelize", indent=False)

num_cores_input = widgets.IntText(
    value=None,  # Default to None, user leaves empty for auto or provides int
    placeholder="Optional: Num Cores",
    description="Num Cores:",
    disabled=True,  # Initially disabled
    style={"description_width": "initial"},
)

core_factor_input = widgets.IntText(
    value=10,
    placeholder="Core Factor",
    description="Core Factor:",
    disabled=True,  # Initially disabled
    style={"description_width": "initial"},
)


# Link parallelize checkbox to enable/disable core inputs
def handle_parallelize_change(change):
    num_cores_input.disabled = not change["new"]
    core_factor_input.disabled = not change["new"]


parallelize_checkbox.observe(handle_parallelize_change, names="value")

# Button
run_button = widgets.Button(
    description="Run Retrosynthesis",
    button_style="success",  # 'success', 'info', 'warning', 'danger' or ''
    tooltip="Click to start the retrosynthesis process",
    icon="cogs",  # (FontAwesome names without the `fa-` prefix)
    layout=widgets.Layout(width="200px"),  # Added layout for width
)

# Output area for run messages
output_area = widgets.Output()

# --- Solution Display Widgets (Defined Statically) ---
# Add Filter checkbox
filter_checkbox = widgets.Checkbox(value=False, description="Filter", indent=False)

fragment_count_input = widgets.IntText(
    value=None,
    placeholder="Enter desired fragment count",
    description="Number of Fragments:",
    style={"description_width": "initial"},
    layout=widgets.Layout(width="300px"),
    disabled=True,  # Initially disabled
)


# Link filter checkbox to enable/disable fragment count input
def handle_filter_change(change):
    fragment_count_input.disabled = not change["new"]
    # Clear the value if disabling the filter
    if not change["new"]:
        fragment_count_input.value = None


filter_checkbox.observe(handle_filter_change, names="value")

display_button = widgets.Button(
    description="Display Solutions",
    button_style="info",
    tooltip="Click to display filtered solutions",
    icon="eye",
    layout=widgets.Layout(width="200px"),
)

# Dropdown to select which solution image to view (initially empty)
solution_dropdown = widgets.Dropdown(
    options=[],
    description="Select Solution:",
    style={"description_width": "initial"},
    value=None,
    disabled=True,  # Start disabled until solutions are loaded
    layout=widgets.Layout(width="auto"),
)

# Output widget to display the selected solution's image
image_display_area = widgets.Output()

# Output area for display messages (errors, status)
solution_output_area = widgets.Output()

# --- SMILES Viewer Widgets ---
fragment_comb_dropdown = widgets.Dropdown(
    options=[],
    description="Select Fragment Comb:",
    style={"description_width": "initial"},
    value=None,
    disabled=True,
    layout=widgets.Layout(width="auto"),
)

# Change from HTML to Output to display images
smiles_display_area = widgets.Output(
    layout=widgets.Layout(min_height="320px", min_width="320px")  # Adjust layout for image
)

smiles_pagination_label = widgets.Label(value="0 of 0")

prev_smiles_button = widgets.Button(
    description="Previous", icon="arrow-left", disabled=True, layout=widgets.Layout(width="100px")
)

next_smiles_button = widgets.Button(
    description="Next", icon="arrow-right", disabled=True, layout=widgets.Layout(width="100px")
)

# --- State ---
# Simple dictionary to hold application state like the solution object
app_state = {
    "retro_solution": None,
    "valid_images_cache": [],  # Cache images for the observer
    "retro_tool": None,  # Store the retro tool instance
    "selected_fragment_comb": None,
    "current_smiles_list": [],
    "current_smiles_index": 0,
    "displayable_solutions": [],  # Store solutions corresponding to valid images
}

# --- Logic ---


def run_retrosynthesis_on_click(b):
    app_state["retro_solution"] = None  # Reset state on new run
    app_state["valid_images_cache"] = []  # Clear image cache
    solution_dropdown.options = []  # Clear dropdown
    solution_dropdown.value = None
    solution_dropdown.disabled = True
    with image_display_area:
        clear_output()
    with solution_output_area:
        clear_output()
    # Reset SMILES viewer state
    fragment_comb_dropdown.options = []
    fragment_comb_dropdown.value = None
    fragment_comb_dropdown.disabled = True
    prev_smiles_button.disabled = True
    next_smiles_button.disabled = True
    smiles_pagination_label.value = "0 of 0"
    with smiles_display_area:
        clear_output()
        # Use logger instead of print
        logger.info("[GUI] Run retrosynthesis to populate fragment combinations.")
    app_state["retro_tool"] = None
    app_state["selected_fragment_comb"] = None
    app_state["current_smiles_list"] = []
    app_state["current_smiles_index"] = 0
    app_state["displayable_solutions"] = []  # Clear displayable solutions

    with output_area:
        clear_output(wait=True)  # Clear previous run output
        # Use logger instead of print
        logger.info("[GUI] Starting retrosynthesis...")

        target = target_smiles_input.value
        fragmenter_name = fragmenter_choice.value
        json_path_str = file_path_input.value
        parallelize = parallelize_checkbox.value
        num_cores = num_cores_input.value  # Will be None if empty
        core_factor = core_factor_input.value

        if not target:
            # Use logger instead of print
            logger.error("[GUI] ERROR: Target SMILES cannot be empty.")
            return

        if not json_path_str:
            # Use logger instead of print
            logger.warning("[GUI] WARNING: Properties JSON path is empty. Proceeding without molecule properties.")
            mol_properties_path = None
        else:
            mol_properties_path = Path(json_path_str)
            if not mol_properties_path.is_file():
                # Use logger instead of print
                logger.error(f"[GUI] ERROR: Properties file not found at {mol_properties_path}")
                # Fallback or further error handling needed
                # For now, setting to None to allow process without properties
                mol_properties_path = None
                # Use logger instead of print
                logger.warning("[GUI] WARNING: Proceeding without molecule properties.")
                # return # Or decide how to proceed

        try:
            # Select fragmenter
            if fragmenter_name == "BRICSFragmenter":
                fragmenter = BRICSFragmenter(target)
            elif fragmenter_name == "rBRICSFragmenter":
                fragmenter = rBRICSFragmenter(target)
            else:
                # Use logger instead of print
                logger.error(f"[GUI] ERROR: Unknown fragmenter type '{fragmenter_name}'")
                return

            # Use logger instead of print
            logger.info(f"[GUI] Using Fragmenter: {fragmenter_name}")
            logger.info(f"[GUI] Using Properties: {mol_properties_path}")
            logger.info(f"[GUI] Parallelize: {parallelize}, Num Cores: {num_cores}, Core Factor: {core_factor}")

            # Initialize Retrosynthesis tool
            # Pass parallelization options
            retro_tool = Retrosynthesis(
                fragmenter,
                original_BBs=None,
                mol_properties_path=mol_properties_path,
                parallelize=parallelize,
                num_cores=num_cores,
                core_factor=core_factor,
            )

            # Run fragmentation and retrosynthesis
            # Use logger instead of print
            logger.info("[GUI] Running fragmentation and retrosynthesis...")
            retro_tool.fragment_retrosynthesis()

            # Process solutions
            retro_solution = RetrosynthesisSolution(retro_tool)
            retro_solution.fill_solutions()

            # Use logger instead of print
            logger.info(f"[GUI] Found {len(retro_solution.solutions)} solution(s).")

            # Store the solution object and retro_tool in the state
            app_state["retro_solution"] = retro_solution
            app_state["retro_tool"] = retro_tool  # Store the tool
            # Use logger instead of print
            logger.info("[GUI] Retrosynthesis complete. Ready to display solutions and browse fragment SMILES.")

            # Populate fragment combination dropdown - REMOVED FROM HERE
            # The dropdown will now be populated when solutions are displayed
            # Reset/disable the dropdown initially
            fragment_comb_dropdown.options = []
            fragment_comb_dropdown.value = None
            fragment_comb_dropdown.disabled = True
            prev_smiles_button.disabled = True
            next_smiles_button.disabled = True
            smiles_pagination_label.value = "0 of 0"
            with smiles_display_area:
                clear_output()
                # Use logger instead of print
                logger.info("[GUI] Display solutions to populate fragment combinations.")

        except Exception as e:
            # Use logger instead of print, include traceback
            logger.error(f"[GUI] ERROR during retrosynthesis: {e}", exc_info=True)


# Attach the function to the button's click event
run_button.on_click(run_retrosynthesis_on_click)


# --- Helper to Update Fragment Combination Dropdown ---
def update_fragment_comb_dropdown(solution):
    """Populates the fragment comb dropdown based on a single solution."""
    # Temporarily unobserve to prevent firing during update
    try:
        fragment_comb_dropdown.unobserve(on_fragment_comb_select, names="value")  # Temporarily unobserve
    except ValueError:
        pass  # Observer might not be attached

    if solution:
        # Extract unique fragment combinations from the *single* selected solution
        # Sort them for consistent order (tuples are sortable)
        unique_combs = sorted(list(set(solution)))
        fragment_comb_dropdown.options = [(str(comb), comb) for comb in unique_combs]
        fragment_comb_dropdown.disabled = False
        with smiles_display_area:
            clear_output()
            # Use logger instead of print
            logger.info("[GUI] Select a fragment combination to view SMILES.")
    else:
        # No solution provided or solution is empty
        fragment_comb_dropdown.options = []
        fragment_comb_dropdown.disabled = True
        with smiles_display_area:
            clear_output()
            # Use logger instead of print
            logger.info("[GUI] No fragment combinations in the selected solution.")

    # Reset SMILES viewer state
    fragment_comb_dropdown.value = None
    prev_smiles_button.disabled = True
    next_smiles_button.disabled = True
    smiles_pagination_label.value = "0 of 0"
    app_state["selected_fragment_comb"] = None
    app_state["current_smiles_list"] = []
    app_state["current_smiles_index"] = 0

    fragment_comb_dropdown.observe(on_fragment_comb_select, names="value")  # Re-register observer


# --- SMILES Viewer Logic ---


def update_smiles_display():
    """Updates the SMILES display area and pagination controls."""
    smiles_list = app_state.get("current_smiles_list", [])
    index = app_state.get("current_smiles_index", 0)
    total_smiles = len(smiles_list)

    with smiles_display_area:
        clear_output(wait=True)  # Clear previous image/message
        if not smiles_list:
            # Use logger instead of print
            logger.info("[GUI] No SMILES found for this combination.")
            smiles_pagination_label.value = "0 of 0"
            prev_smiles_button.disabled = True
            next_smiles_button.disabled = True
            return

        current_smiles = smiles_list[index]

        try:
            mol = Chem.MolFromSmiles(current_smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))  # Generate PIL image
                # Convert PIL image to PNG bytes
                bio = io.BytesIO()
                img.save(bio, format="PNG")
                image_widget = widgets.Image(value=bio.getvalue(), format="png", width=300, height=300)
                display(image_widget)
            else:
                # Use logger instead of print
                logger.warning(f"[GUI] Invalid SMILES: {current_smiles}")
        except Exception as e:
            # Use logger instead of print
            logger.error(f"[GUI] Error generating image for SMILES {current_smiles}: {e}")

    # Update pagination label and button states outside the 'with' block
    # as they are separate widgets
    smiles_pagination_label.value = f"{index + 1} of {total_smiles}"
    prev_smiles_button.disabled = index == 0
    next_smiles_button.disabled = index >= total_smiles - 1


def on_fragment_comb_select(change):
    """Handles selection changes in the fragment combination dropdown."""
    if change["type"] == "change" and change["name"] == "value":
        selected_comb = change["new"]
        retro_tool = app_state.get("retro_tool")

        if selected_comb is None or retro_tool is None or not hasattr(retro_tool, "comb_bbs_dict"):
            app_state["selected_fragment_comb"] = None
            app_state["current_smiles_list"] = []
            app_state["current_smiles_index"] = 0
        else:
            # Ensure smiles_list is a list for indexing
            smiles_data = retro_tool.comb_bbs_dict.get(selected_comb, set())  # Default to empty set
            smiles_list = sort_by_heavy_atoms(list(smiles_data))  # Convert set (or other iterable) to list
            app_state["selected_fragment_comb"] = selected_comb
            app_state["current_smiles_list"] = smiles_list
            app_state["current_smiles_index"] = 0  # Reset index on new selection

        update_smiles_display()  # Update the display


def on_prev_smiles_click(b):
    """Handles clicks on the 'Previous' SMILES button."""
    if app_state["current_smiles_index"] > 0:
        app_state["current_smiles_index"] -= 1
        update_smiles_display()


def on_next_smiles_click(b):
    """Handles clicks on the 'Next' SMILES button."""
    smiles_list = app_state.get("current_smiles_list", [])
    if app_state["current_smiles_index"] < len(smiles_list) - 1:
        app_state["current_smiles_index"] += 1
        update_smiles_display()


# Attach handlers
fragment_comb_dropdown.observe(on_fragment_comb_select, names="value")
prev_smiles_button.on_click(on_prev_smiles_click)
next_smiles_button.on_click(on_next_smiles_click)


# --- Solution Display Logic ---
def display_solutions_on_click(b):
    # Clear previous display messages and image
    with solution_output_area:
        clear_output(wait=True)
        # Use logger instead of print
        logger.info("[GUI] Attempting to display solutions...")

    with image_display_area:
        clear_output(wait=True)

    # Reset solution dropdown state
    try:
        solution_dropdown.unobserve(on_solution_select, names="value")  # Temporarily unobserve
    except ValueError:
        pass  # Observer likely wasn't attached yet
    solution_dropdown.options = []
    solution_dropdown.value = None
    solution_dropdown.disabled = True
    app_state["valid_images_cache"] = []  # Clear cache
    app_state["displayable_solutions"] = []  # Clear displayable solutions

    # --- Reset SMILES viewer state ---
    try:
        fragment_comb_dropdown.unobserve(on_fragment_comb_select, names="value")  # Temporarily unobserve
    except ValueError:
        pass  # Observer likely wasn't attached yet
    fragment_comb_dropdown.options = []
    fragment_comb_dropdown.value = None
    fragment_comb_dropdown.disabled = True
    prev_smiles_button.disabled = True
    next_smiles_button.disabled = True
    smiles_pagination_label.value = "0 of 0"
    with smiles_display_area:
        clear_output()
        # Use logger instead of print
        logger.info("[GUI] Select a solution to view its fragment combinations.")
    app_state["selected_fragment_comb"] = None
    app_state["current_smiles_list"] = []
    app_state["current_smiles_index"] = 0
    # --- End Reset SMILES viewer state ---

    retro_solution = app_state.get("retro_solution")
    use_filter = filter_checkbox.value
    fragment_count = fragment_count_input.value if use_filter else None  # Get value only if filter is checked

    if retro_solution is None:
        msg = "No retrosynthesis results available. Please run retrosynthesis first."
        logger.warning(msg)
        with solution_output_area:
            # Use logger instead of print
            logger.warning(f"[GUI] {msg}")
        return

    # Apply filter only if checkbox is checked and value is valid
    if use_filter and (fragment_count is None or not isinstance(fragment_count, int) or fragment_count <= 0):
        msg = "Filter checkbox is checked, but fragment count is invalid. Displaying all solutions."
        logger.warning(f"[GUI] {msg}")
        with solution_output_area:
            # Use logger instead of print
            logger.warning(f"[GUI] {msg}")
        filtered_solutions = retro_solution.solutions
    elif use_filter:
        # Checkbox is checked and value is valid - apply filter
        filtered_solutions = [sol for sol in retro_solution.solutions if len(sol) == fragment_count]
        msg = f"Filtering for solutions with exactly {fragment_count} fragments."
        logger.info(f"[GUI] {msg}")
        with solution_output_area:
            # Use logger instead of print
            logger.info(f"[GUI] {msg}")
    else:
        # Filter checkbox is unchecked - display all solutions
        msg = "Filter checkbox is unchecked. Displaying all available solutions."
        # Use logger instead of print
        logger.info(f"[GUI] {msg}")
        with solution_output_area:
            # Use logger instead of print
            logger.info(f"[GUI] {msg}")
        filtered_solutions = retro_solution.solutions

    if not filtered_solutions:
        msg = f"No solutions found matching the criteria (count: {fragment_count})."
        logger.info(f"[GUI] {msg}")
        with solution_output_area:
            # Use logger instead of print
            logger.info(f"[GUI] {msg}")
        return

    # Use logger instead of print
    logger.info(f"[GUI] Visualizing {len(filtered_solutions)} solution(s)...")
    with solution_output_area:
        # Use logger instead of print
        logger.info(f"[GUI] Visualizing {len(filtered_solutions)} solution(s)...")

    try:
        # Generate one image per solution in filtered_solutions
        solution_images = retro_solution.visualize_solutions(filtered_solutions)

        if not solution_images:
            msg = "Visualization did not produce any images."
            logger.info(f"[GUI] {msg}")
            with solution_output_area:
                # Use logger instead of print
                logger.info(f"[GUI] {msg}")
            return

        # Filter out None values and store valid images and corresponding solutions
        valid_images = []
        displayable_solutions = []  # Temporary list for this run
        original_indices = []  # Keep track of original index for labeling
        for i, img in enumerate(solution_images):
            if img:
                valid_images.append(img)
                # Store the actual solution corresponding to this valid image
                displayable_solutions.append(filtered_solutions[i])
                # Assuming order is preserved to find original index
                original_solution_index = retro_solution.solutions.index(filtered_solutions[i])
                original_indices.append(original_solution_index)
            else:
                # Use logger instead of print
                logger.warning(f"[GUI] Null image returned by visualize_solutions for filtered solution index {i}")

        num_valid_images = len(valid_images)
        app_state["valid_images_cache"] = valid_images  # Store images for the observer
        app_state["displayable_solutions"] = displayable_solutions  # Store solutions for the observer

        if num_valid_images == 0:
            msg = "Visualization did not produce any valid images."
            logger.info(f"[GUI] {msg}")
            with solution_output_area:
                # Use logger instead of print
                logger.info(f"[GUI] {msg}")
            return

        # Use logger instead of print
        logger.info(f"[GUI] Generated {num_valid_images} image(s) for interactive display.")
        with solution_output_area:
            # Use logger instead of print
            logger.info(f"[GUI] Generated {num_valid_images} image(s). Use dropdown to view.")

        # --- Update Interactive Display Widgets ---

        # Update Solution Dropdown options
        # The value `i` will now be the index into valid_images and displayable_solutions
        dropdown_options = [(f"Solution {original_indices[i]+1}", i) for i in range(num_valid_images)]
        solution_dropdown.options = dropdown_options
        solution_dropdown.value = 0  # Default to showing the first valid solution
        solution_dropdown.disabled = False
        solution_dropdown.observe(on_solution_select, names="value")  # Re-register observer

        # --- Populate Fragment Combination Dropdown for the INITIALLY selected solution ---
        if num_valid_images > 0:
            initial_solution = app_state["displayable_solutions"][0]
            update_fragment_comb_dropdown(initial_solution)
        else:
            # No valid solutions, ensure dropdown is empty/disabled
            update_fragment_comb_dropdown(None)

        # --- End Initial Fragment Combination Population ---

        # Display the initial image
        with image_display_area:
            clear_output(wait=True)
            display(valid_images[0])  # Display the first valid image initially

        # Use logger instead of print
        logger.info("[GUI] Solutions displayed.")

    except Exception as e:
        # Use logger instead of print
        logger.error(f"[GUI] An error occurred during solution visualization: {e}", exc_info=True)
        with solution_output_area:
            # Use logger instead of print
            logger.error(f"[GUI] An error occurred during visualization: {e}")


# Attach the function to the display button's click event
display_button.on_click(display_solutions_on_click)

# --- Layout ---

# Combine widgets into a VBox (vertical box layout)
# Group parallelization options
parallel_options = widgets.VBox([parallelize_checkbox, num_cores_input, core_factor_input])

# VBox for the solution display section
solution_display_section = widgets.VBox(
    [
        widgets.HTML("<hr><h4>Display Solutions:</h4>"),  # Separator and header
        filter_checkbox,
        fragment_count_input,
        display_button,
        solution_output_area,  # Area for messages related to display
        solution_dropdown,  # Static dropdown
        image_display_area,  # Static image output area
    ]
)

# VBox for the SMILES viewer section
smiles_pagination_controls = widgets.HBox([prev_smiles_button, smiles_pagination_label, next_smiles_button])
smiles_viewer_section = widgets.VBox(
    [
        widgets.HTML("<hr><h4>Browse Fragment Combination SMILES:</h4>"),
        fragment_comb_dropdown,
        smiles_display_area,
        smiles_pagination_controls,
    ]
)

gui_layout = widgets.VBox(
    [
        widgets.HTML("<h2>FragmentRetro Retrosynthesis GUI</h2>"),
        target_smiles_input,
        fragmenter_choice,
        file_path_input,
        widgets.HTML("<h4>Parallelization Options:</h4>"),  # Add a header for clarity
        parallel_options,  # Add the group here
        run_button,
        output_area,  # Area for run messages
        solution_display_section,  # Add the whole solution display section here
        smiles_viewer_section,  # Add the SMILES viewer section
    ]
)


# Function to display the GUI in a notebook
def display_gui(smiles: str | None = None):
    """Displays the Retrosynthesis GUI. Optionally pre-fills the target SMILES."""
    if smiles:
        target_smiles_input.value = smiles
    # Ensure initial state is clean
    solution_dropdown.options = []
    solution_dropdown.value = None
    solution_dropdown.disabled = True
    with image_display_area:
        clear_output()
    with solution_output_area:
        clear_output()
    with output_area:
        clear_output()
    # Reset SMILES viewer state
    fragment_comb_dropdown.options = []
    fragment_comb_dropdown.value = None
    fragment_comb_dropdown.disabled = True
    prev_smiles_button.disabled = True
    next_smiles_button.disabled = True
    smiles_pagination_label.value = "0 of 0"
    with smiles_display_area:
        clear_output()
        # Use logger instead of print
        logger.info("[GUI] Run retrosynthesis to view fragment combination SMILES.")
    app_state["retro_tool"] = None
    app_state["selected_fragment_comb"] = None
    app_state["current_smiles_list"] = []
    app_state["current_smiles_index"] = 0
    app_state["displayable_solutions"] = []  # Clear displayable solutions

    display(gui_layout)


# --- Observer for Dropdown --- (Defined globally once)
def on_solution_select(change):
    if change["type"] == "change" and change["name"] == "value":
        selected_index = change[
            "new"
        ]  # The 'new' key holds the selected value (index in valid_images_cache/displayable_solutions)
        valid_images_cache = app_state.get("valid_images_cache", [])
        displayable_solutions = app_state.get("displayable_solutions", [])

        # Display the selected image
        if selected_index is not None and 0 <= selected_index < len(valid_images_cache):
            with image_display_area:
                clear_output(wait=True)  # Clear previous image
                display(valid_images_cache[selected_index])  # Directly display the selected image object

            # Update the fragment comb dropdown based on the selected solution
            if selected_index < len(displayable_solutions):
                selected_solution = displayable_solutions[selected_index]
                update_fragment_comb_dropdown(selected_solution)
            else:
                # Should not happen if lists are kept in sync, but handle defensively
                update_fragment_comb_dropdown(None)

        elif selected_index is None:
            with image_display_area:
                clear_output(wait=True)
                # Use logger instead of print
                logger.info("[GUI] No solution selected or available.")
            # Clear the fragment comb dropdown as well
            update_fragment_comb_dropdown(None)
