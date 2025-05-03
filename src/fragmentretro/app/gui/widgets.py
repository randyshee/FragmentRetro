import ipywidgets as widgets

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

# --- Solution Display Widgets ---
# Add Filter checkbox
filter_checkbox = widgets.Checkbox(value=False, description="Filter Solutions (by fragment count):", indent=False)

fragment_count_input = widgets.IntText(
    value=None,
    placeholder="Enter desired fragment count",
    description="Number of Fragments:",
    style={"description_width": "initial"},
    layout=widgets.Layout(width="300px"),
    disabled=True,  # Initially disabled
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

# NEW: Button to sort SMILES
sort_smiles_button = widgets.Button(
    description="Sort by Heavy Atoms",
    icon="sort",
    tooltip="Click to sort the displayed SMILES list by heavy atom count",
    disabled=True,
    layout=widgets.Layout(width="180px"),
)
