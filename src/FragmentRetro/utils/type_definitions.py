"""Type definitions for FragmentRetro package."""

from typing import TypedDict

from typing_extensions import TypeAlias

# Bond type (from BRICS): ((atom1_idx, atom2_idx), (type1, type2))
# Example: ((3, 2), ('3', '4'))
BondType: TypeAlias = tuple[tuple[int, int], tuple[str, str]]

# Atom mapping type: tuple of atom indices
# Example: (0, 1, 2, 7)
AtomMappingType: TypeAlias = tuple[int, ...]

# Combination type: tuple of fragment indices (needs to be a tuple for hashing)
CombType: TypeAlias = tuple[int]

# Solution type: list of combination types
SolutionType: TypeAlias = list[CombType]

# Building block type: set of strings (SMILES)
BBsType: TypeAlias = set[str]

# Store the combinations for each stage
# Example: {1: [(0, 1), (2, 3)]} ({stage number: list[CombType]})
# This can hold all combinations or effective combinations
StageCombDictType: TypeAlias = dict[int, list[CombType]]

# Store the building blocks for each combination
CombBBsDictType: TypeAlias = dict[CombType, BBsType]


class MolProperties(TypedDict, total=True):
    """Type definition for molecular properties."""

    cano_smiles: str
    num_heavy_atoms: int
    num_rings: int
    pfp: list[int]  # Pattern Fingerprint
