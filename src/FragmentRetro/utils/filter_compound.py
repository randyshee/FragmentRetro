import json
from pathlib import Path
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm

from FragmentRetro.utils.helpers import canonicalize_smiles
from FragmentRetro.utils.logging_config import logger
from FragmentRetro.utils.type_definitions import (
    BBsType,
    FilterIndicesType,
    MolProperties,
)


def get_mol_properties(smiles: str, fpSize: int = 2048) -> MolProperties:
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

    pfp = list(Chem.rdmolops.PatternFingerprint(mol, fpSize=fpSize).GetOnBits())

    return {
        "cano_smiles": cano_smiles,
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "pfp": pfp,
    }


def precompute_properties(smiles_list: list[str], output_path: Path, fpSize: int = 2048) -> None:
    """Calculates molecular properties for a list of SMILES strings and saves them to a JSON file.

    Args:
        smiles_list: A list of SMILES strings.
        output_path: The path to the output JSON file.
    """
    results = []
    for smiles in tqdm(smiles_list, desc="Precomputing molecular properties"):
        try:
            mol_properties = get_mol_properties(smiles, fpSize=fpSize)
            results.append(mol_properties)
        except ValueError as e:
            logger.error(f"Error processing SMILES '{smiles}': {e} during precompute_properties")
            continue

    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)


class CompoundFilter:
    """A class for filtering compounds based on precomputed molecular properties."""

    def __init__(self, mol_properties_path: Path, fpSize: int = 2048):
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
        self.fpSize = fpSize

        self._load_mol_properties()
        self._create_numpy_arrays()

    def _load_mol_properties(self) -> None:
        """Loads molecular properties from the JSON file."""
        with open(self.mol_properties_path, "r") as f:
            mol_properties_list = json.load(f)

        self.len_BBs = len(mol_properties_list)

        for mol_props in mol_properties_list:
            self.cano_smiles_list.append(mol_props["cano_smiles"])
            self.num_heavy_atoms_list.append(mol_props["num_heavy_atoms"])
            self.num_rings_list.append(mol_props["num_rings"])
            self.pfp_len_list.append(len(mol_props["pfp"]))
            self.pfp_list.append(mol_props["pfp"])

    def _create_numpy_arrays(self) -> None:
        """Creates NumPy arrays for faster filtering."""
        self.num_heavy_atoms_array = np.array(self.num_heavy_atoms_list)
        self.num_rings_array = np.array(self.num_rings_list)
        self.pfp_len_array = np.array(self.pfp_len_list)

        # Create a boolean NumPy array for PFP bits
        self.pfp_bit_array = np.zeros((len(self.pfp_list), self.fpSize), dtype=bool)
        for i, pfp in enumerate(self.pfp_list):
            self.pfp_bit_array[i, pfp] = True

    def filter_compounds(
        self, smiles: str, prefiltered_indices: Optional[FilterIndicesType] = None
    ) -> FilterIndicesType:
        """Filters compounds based on a query SMILES string and prefiltered indices.

        Args:
            smiles: The query SMILES string.
            prefiltered_indices: A list of prefiltered indices.

        Returns:
            A list of indices of the compounds that pass the filter.
        """
        try:
            mol_properties = get_mol_properties(smiles, fpSize=self.fpSize)
        except ValueError as e:
            print(f"Invalid SMILES: {e}")
            return []

        logger.info(f"Filtering BBs for {smiles}")

        num_heavy_atoms = mol_properties["num_heavy_atoms"]
        num_rings = mol_properties["num_rings"]
        pfp = mol_properties["pfp"]
        pfp_len = len(pfp)

        query_pfp_bit_array = np.zeros(self.fpSize, dtype=bool)
        query_pfp_bit_array[pfp] = True

        # Filtering based on molecular properties
        indices_array = np.where(
            (self.num_heavy_atoms_array >= num_heavy_atoms)
            & (self.num_rings_array >= num_rings)
            & (self.pfp_len_array >= pfp_len)
        )[0]
        if prefiltered_indices is not None:
            indices_array = np.intersect1d(indices_array, prefiltered_indices)

        # check pfp of query is a subset of pfp of filtered compounds
        if indices_array.size == 0:
            filtered_indices = []
        else:
            filtered_indices = indices_array[
                np.all(self.pfp_bit_array[indices_array][:, query_pfp_bit_array], axis=1)
            ].tolist()

        if prefiltered_indices is None:
            logger.info(f"Originally {self.len_BBs} BBs, filtered down to {len(filtered_indices)}")
        else:
            logger.info(f"Originally {len(prefiltered_indices)} BBs, filtered down to {len(filtered_indices)}")

        return filtered_indices

    def get_filtered_BBs(
        self, smiles: str, prefiltered_indices: Optional[FilterIndicesType] = None
    ) -> tuple[FilterIndicesType, BBsType]:
        """Filters building blocks based on a query SMILES string and prefiltered indices.

        This method filters the building blocks based on the properties
        of the provided SMILES string and a list of prefiltered indices.
        It uses the `filter_compounds` method to get a list of indices
        that pass the filter, and then returns a set of the corresponding
        canonical SMILES strings.

        Args:
            smiles: The query SMILES string.
            prefiltered_indices: A list of prefiltered indices.

        Returns:
            tuple[list[int], BBsType]: A tuple containing a list of indices
            of the compounds that pass the filter and a set of canonical
            SMILES strings of the building blocks that pass the filter.
        """
        filtered_indices = self.filter_compounds(smiles, prefiltered_indices)
        return filtered_indices, set(self.cano_smiles_list[i] for i in filtered_indices)
