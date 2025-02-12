import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from FragmentRetro.utils.helpers import canonicalize_smiles
from FragmentRetro.utils.type_definitions import MolProperties


def get_mol_properties(smiles: str) -> MolProperties:
    """Given a SMILES string, returns a dictionary containing molecular properties.

    Args:
        smiles: The SMILES string.

    Returns:
        A dictionary containing the following keys:
            - 'num_heavy_atoms': Number of heavy atoms
            - 'num_rings': Number of rings
            - 'pfp': Pattern fingerprint bits

    Raises:
        ValueError: If the SMILES string is invalid and cannot be converted
            to an RDKit molecule.
    """
    cano_smiles = canonicalize_smiles(smiles)
    mol = Chem.MolFromSmiles(cano_smiles)
    # solve C++ signature problems?
    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)

    pfp = list(Chem.rdmolops.PatternFingerprint(mol).GetOnBits())

    return {
        "cano_smiles": cano_smiles,
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "pfp": pfp,
    }


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


class CompoundFilter:
    """A class for filtering compounds based on precomputed molecular properties."""

    def __init__(self, mol_properties_path: Path):
        """
        Initializes the CompoundFilter with molecular properties loaded from a JSON file.

        Args:
            mol_properties_path: Path to the JSON file containing molecular properties.
        """
        self.mol_properties_path = mol_properties_path
        self.cano_smiles_list: list[str] = []
        self.num_heavy_atoms_list: list[int] = []
        self.num_rings_list: list[int] = []
        self.pfp_len_list: list[int] = []
        self.pfp_list: list[list[int]] = []

        self._load_mol_properties()
        self._create_numpy_arrays()

    def _load_mol_properties(self):
        """Loads molecular properties from the JSON file."""
        with open(self.mol_properties_path, "r") as f:
            mol_properties_list = json.load(f)

        for mol_props in mol_properties_list:
            self.cano_smiles_list.append(mol_props["cano_smiles"])
            self.num_heavy_atoms_list.append(mol_props["num_heavy_atoms"])
            self.num_rings_list.append(mol_props["num_rings"])
            self.pfp_len_list.append(len(mol_props["pfp"]))
            self.pfp_list.append(mol_props["pfp"])

    def _create_numpy_arrays(self):
        """Creates NumPy arrays for faster filtering."""
        self.num_heavy_atoms_array = np.array(self.num_heavy_atoms_list)
        self.num_rings_array = np.array(self.num_rings_list)
        self.pfp_len_array = np.array(self.pfp_len_list)

        # Create a set of all unique PFP bits
        all_pfp_bits = set()
        for pfp in self.pfp_list:
            all_pfp_bits.update(pfp)
        self.all_pfp_bits = sorted(list(all_pfp_bits))

        # Create a boolean NumPy array for PFP bits
        self.pfp_bit_array = np.zeros((len(self.pfp_list), len(self.all_pfp_bits)), dtype=bool)
        for i, pfp in enumerate(self.pfp_list):
            self.pfp_bit_array[i, [self.all_pfp_bits.index(bit) for bit in pfp]] = True

    def filter_compounds(self, smiles: str) -> list[int]:
        """Filters compounds based on a query SMILES string.

        Args:
            smiles: The query SMILES string.

        Returns:
            A list of indices of the compounds that pass the filter.
        """
        try:
            mol_properties = get_mol_properties(smiles)
        except ValueError as e:
            print(f"Invalid SMILES: {e}")
            return []

        num_heavy_atoms = mol_properties["num_heavy_atoms"]
        num_rings = mol_properties["num_rings"]
        pfp = mol_properties["pfp"]
        pfp_len = len(pfp)

        # Filtering based on molecular properties
        indices = np.where(
            (self.num_heavy_atoms_array >= num_heavy_atoms)
            & (self.num_rings_array >= num_rings)
            & (self.pfp_len_array >= pfp_len)
        )[0].tolist()

        # # Filter based on PFP bits
        # query_pfp_bits = set(pfp)
        # filtered_indices = []
        # for i in indices:
        #     stock_pfp_bits = set(self.pfp_list[i])
        #     if query_pfp_bits.issubset(stock_pfp_bits):
        #         filtered_indices.append(i)

        # check pfpbit matches
        query_indices = [self.all_pfp_bits.index(bit) for bit in pfp if bit in self.all_pfp_bits]
        has_bits_matched = np.all(self.pfp_bit_array[indices][:, query_indices], axis=1)
        filtered_indices = np.array(indices)[has_bits_matched].tolist()

        return filtered_indices
