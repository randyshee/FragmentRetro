import re
from typing import cast

from rdkit import Chem


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


def replace_dummy_atoms_regex(smiles: str) -> str:
    """Replaces dummy atoms ('*') in a SMILES string with explicit hydrogen ('H') using regex.

    Args:
        smiles: The SMILES string containing dummy atoms.

    Returns:
        A SMILES string where '[{int}*]' is replaced with '[H]'.
    """
    return canonicalize_smiles(re.sub(r"\[\d*\*\]", "[H]", smiles))


def remove_indices_before_dummy(smiles: str) -> str:
    """Removes indices before asterisks (*) in a SMILES string using regex.

    Args:
        smiles: The SMILES string containing indices before asterisks.

    Returns:
        The SMILES string with indices before asterisks removed.
    """
    return re.sub(r"\[\d*\*\]", "[*]", smiles)
