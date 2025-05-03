import ipywidgets as widgets

from app.gui.widgets import (
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
        sort_smiles_button,  # Add sort button to layout
    ]
)

gui_layout = widgets.VBox(
    [
        widgets.HTML("<h2>FragmentRetro GUI</h2>"),
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
