import json
from pathlib import Path

from FragmentRetro.utils.helpers import get_mol_properties


def precompute_properties(smiles_list: list[str], output_path: Path) -> None:
    """Calculates molecular properties for a list of SMILES strings and saves them to a JSON file.

    Args:
        smiles_list: A list of SMILES strings.
        output_path: The path to the output JSON file.
    """
    results = []
    for smiles in smiles_list:
        try:
            mol_properties = get_mol_properties(smiles)
            results.append(mol_properties)
        except ValueError as e:
            print(f"Error processing SMILES '{smiles}': {e}")
            continue

    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)
