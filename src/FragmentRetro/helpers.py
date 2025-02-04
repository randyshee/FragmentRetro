from rdkit import Chem


def count_heavy_atoms(smiles: str) -> int:
    """Counts the number of heavy atoms in a SMILES string.

    Heavy atoms are defined as atoms with an atomic number greater than 1,
    excluding hydrogen.

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
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


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
