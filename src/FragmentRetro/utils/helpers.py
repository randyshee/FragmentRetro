from typing import cast

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from FragmentRetro.utils.type_definitions import MolProperties


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalizes a SMILES string using RDKit.

    Args:
        smiles: The SMILES string to canonicalize.

    Returns:
        The canonicalized SMILES string.

    Raises:
        ValueError: If the SMILES string cannot be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")
    return cast(str, Chem.MolToSmiles(mol))


def count_heavy_atoms(smiles: str) -> int:
    """Counts the number of heavy atoms in a SMILES string.

    Args:
        smiles: The SMILES string representing the chemical structure.

    Returns:
        The number of heavy atoms in the molecule.

    Raises:
        ValueError: If the SMILES string is invalid and cannot be converted
            to an RDKit molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return cast(int, mol.GetNumHeavyAtoms())


def sort_by_heavy_atoms(smiles_list: list[str]) -> list[str]:
    """Sorts a list of SMILES strings by the number of heavy atoms.

    The sorting is done in ascending order, i.e., the molecule with the fewest
    heavy atoms will come first.

    Args:
        smiles_list: The list of SMILES strings to be sorted.

    Returns:
        The sorted list of SMILES strings.
    """
    return sorted(smiles_list, key=count_heavy_atoms)


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
