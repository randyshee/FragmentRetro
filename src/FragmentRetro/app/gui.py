from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output, display

from FragmentRetro.fragmenter import BRICSFragmenter, rBRICSFragmenter
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution
from FragmentRetro.utils.logging_config import logger

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
fragment_count_input = widgets.IntText(
    value=None,
    placeholder="Enter desired fragment count",
    description="Number of Fragments:",
    style={"description_width": "initial"},
    layout=widgets.Layout(width="300px"),
)

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

# --- State ---
# Simple dictionary to hold application state like the solution object
app_state = {
    "retro_solution": None,
    "valid_images_cache": [],  # Cache images for the observer
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

    with output_area:
        clear_output(wait=True)  # Clear previous run output
        logger.info("Starting retrosynthesis...")
        print("Starting retrosynthesis...")  # Also print to output area

        target = target_smiles_input.value
        fragmenter_name = fragmenter_choice.value
        json_path_str = file_path_input.value
        parallelize = parallelize_checkbox.value
        num_cores = num_cores_input.value  # Will be None if empty
        core_factor = core_factor_input.value

        if not target:
            logger.error("Target SMILES cannot be empty.")
            print("ERROR: Target SMILES cannot be empty.")
            return

        if not json_path_str:
            logger.warning("Properties JSON path is empty. Proceeding without molecule properties.")
            print("WARNING: Properties JSON path is empty. Proceeding without molecule properties.")
            mol_properties_path = None
        else:
            mol_properties_path = Path(json_path_str)
            if not mol_properties_path.is_file():
                logger.error(f"Properties file not found at {mol_properties_path}")
                print(f"ERROR: Properties file not found at {mol_properties_path}")
                # Fallback or further error handling needed
                # For now, setting to None to allow process without properties
                mol_properties_path = None
                logger.warning("Proceeding without molecule properties.")
                print("WARNING: Proceeding without molecule properties.")
                # return # Or decide how to proceed

        try:
            # Select fragmenter
            if fragmenter_name == "BRICSFragmenter":
                fragmenter = BRICSFragmenter(target)
            elif fragmenter_name == "rBRICSFragmenter":
                fragmenter = rBRICSFragmenter(target)
            else:
                logger.error(f"Unknown fragmenter type '{fragmenter_name}'")
                print(f"ERROR: Unknown fragmenter type '{fragmenter_name}'")
                return

            logger.info(f"Using Fragmenter: {fragmenter_name}")
            print(f"Using Fragmenter: {fragmenter_name}")
            logger.info(f"Using Properties: {mol_properties_path}")
            print(f"Using Properties: {mol_properties_path}")
            logger.info(f"Parallelize: {parallelize}, Num Cores: {num_cores}, Core Factor: {core_factor}")
            print(f"Parallelize: {parallelize}, Num Cores: {num_cores}, Core Factor: {core_factor}")

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
            logger.info("Running fragmentation and retrosynthesis...")
            print("Running fragmentation and retrosynthesis...")
            retro_tool.fragment_retrosynthesis()

            # Process solutions
            retro_solution = RetrosynthesisSolution(retro_tool)
            retro_solution.fill_solutions()

            logger.info(f"Found {len(retro_solution.solutions)} solution(s).")
            print(f"Found {len(retro_solution.solutions)} solution(s).")

            # Store the solution object in the state
            app_state["retro_solution"] = retro_solution
            print("Retrosynthesis complete. Ready to display solutions.")

        except Exception as e:
            logger.error(f"An error occurred during retrosynthesis: {e}", exc_info=True)  # Include traceback
            print(f"ERROR during retrosynthesis: {e}")


# Attach the function to the button's click event
run_button.on_click(run_retrosynthesis_on_click)


# --- Observer for Dropdown --- (Defined globally once)
def on_solution_select(change):
    if change["type"] == "change" and change["name"] == "value":
        selected_index = change["new"]  # The 'new' key holds the selected value (index in valid_images_cache)
        # Ensure the index is valid before updating the image
        valid_images_cache = app_state.get("valid_images_cache", [])
        if selected_index is not None and 0 <= selected_index < len(valid_images_cache):
            with image_display_area:
                clear_output(wait=True)  # Clear previous image
                display(valid_images_cache[selected_index])  # Directly display the selected image object
        elif selected_index is None:
            with image_display_area:
                clear_output(wait=True)
                print("No solution selected or available.")


# Register the observer function with the dropdown widget
solution_dropdown.observe(on_solution_select, names="value")


# --- Solution Display Logic ---
def display_solutions_on_click(b):
    # Clear previous display messages and image
    with solution_output_area:
        clear_output(wait=True)
        logger.info("Attempting to display solutions...")
        print("Attempting to display solutions...")

    with image_display_area:
        clear_output(wait=True)

    # Reset dropdown state
    solution_dropdown.unobserve(on_solution_select, names="value")  # Temporarily unobserve
    solution_dropdown.options = []
    solution_dropdown.value = None
    solution_dropdown.disabled = True
    app_state["valid_images_cache"] = []  # Clear cache

    retro_solution = app_state.get("retro_solution")
    fragment_count = fragment_count_input.value

    if retro_solution is None:
        msg = "No retrosynthesis results available. Please run retrosynthesis first."
        logger.warning(msg)
        with solution_output_area:
            print(msg)
        return

    if fragment_count is None or not isinstance(fragment_count, int) or fragment_count <= 0:
        logger.warning("Invalid fragment count. Displaying all solutions.")
        with solution_output_area:
            print("Invalid fragment count input. Displaying all available solutions.")
        filtered_solutions = retro_solution.solutions
    else:
        filtered_solutions = [sol for sol in retro_solution.solutions if len(sol) == fragment_count]
        msg = f"Filtering for solutions with exactly {fragment_count} fragments."
        logger.info(msg)
        with solution_output_area:
            print(msg)

    if not filtered_solutions:
        msg = f"No solutions found matching the criteria (count: {fragment_count})."
        logger.info(msg)
        with solution_output_area:
            print(msg)
        return

    logger.info(f"Visualizing {len(filtered_solutions)} solution(s)...")
    with solution_output_area:
        print(f"Visualizing {len(filtered_solutions)} solution(s)...")

    try:
        # Generate one image per solution in filtered_solutions
        solution_images = retro_solution.visualize_solutions(filtered_solutions)

        if not solution_images:
            msg = "Visualization did not produce any images."
            logger.info(msg)
            with solution_output_area:
                print(msg)
            return

        # Filter out None values and store valid images with their original indices
        valid_images = []
        original_indices = []  # Keep track of original index for labeling
        for i, img in enumerate(solution_images):
            if img:
                valid_images.append(img)
                # Find original index from filtered_solutions based on content or assume order preservation
                # This assumes visualize_solutions preserves the order of filtered_solutions
                # To be robust, we might need to match based on solution content if order isn't guaranteed
                # For now, assuming order is preserved:
                original_solution_index = retro_solution.solutions.index(filtered_solutions[i])
                original_indices.append(original_solution_index)
            else:
                logger.warning(f"Null image returned by visualize_solutions for filtered solution index {i}")

        num_valid_images = len(valid_images)
        app_state["valid_images_cache"] = valid_images  # Store images for the observer

        if num_valid_images == 0:
            msg = "Visualization did not produce any valid images."
            logger.info(msg)
            with solution_output_area:
                print(msg)
            return

        logger.info(f"Generated {num_valid_images} image(s) for interactive display.")
        with solution_output_area:
            print(f"Generated {num_valid_images} image(s). Use dropdown to view.")

        # --- Update Interactive Display Widgets ---

        # Update Dropdown options
        dropdown_options = [(f"Solution {original_indices[i]+1}", i) for i in range(num_valid_images)]
        solution_dropdown.options = dropdown_options
        solution_dropdown.value = 0  # Default to showing the first valid solution
        solution_dropdown.disabled = False

        # Re-register the observer
        solution_dropdown.observe(on_solution_select, names="value")

        # Display the initial image
        with image_display_area:
            clear_output(wait=True)
            display(valid_images[0])  # Display the first valid image initially

        logger.info("Solutions displayed.")

    except Exception as e:
        logger.error(f"An error occurred during solution visualization: {e}", exc_info=True)
        with solution_output_area:
            print(f"An error occurred during visualization: {e}")


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
        fragment_count_input,
        display_button,
        solution_output_area,  # Area for messages related to display
        solution_dropdown,  # Static dropdown
        image_display_area,  # Static image output area
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
    display(gui_layout)
