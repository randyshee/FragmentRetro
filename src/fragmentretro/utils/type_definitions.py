"""Type definitions for FragmentRetro package."""

from typing import TypedDict

# Bond type (from BRICS): ((atom1_idx, atom2_idx), (type1, type2))
# Example: ((3, 2), ('3', '4'))
type BondType = tuple[tuple[int, int], tuple[str, str]]

# Atom mapping type: tuple of atom indices
# Example: (0, 1, 2, 7)
type AtomMappingType = tuple[int, ...]

# Combination type: tuple of fragment indices (needs to be a tuple for hashing)
type CombType = tuple[int]

# Solution type: list of combination types
type SolutionType = list[CombType]

# Building block type: set of strings (SMILES)
type BBsType = set[str]

# Store the combinations for each stage
# Example: {1: [(0, 1), (2, 3)]} ({stage number: list[CombType]})
# This can hold all combinations or effective combinations
type StageCombDictType = dict[int, list[CombType]]

# Store the building blocks for each combination
type CombBBsDictType = dict[CombType, BBsType]

# Store the combination indices and building blocks for each fragment SMILES (can have dummy atoms)
type FragmentBBsDictType = dict[str, tuple[CombType, BBsType]]

# Store the filtered indices
type FilterIndicesType = list[int]

# Store the filtered indices (from CompoundFilter) for each combination
type CombFilterIndicesDictType = dict[CombType, FilterIndicesType]


class MolProperties(TypedDict, total=True):
    """Type definition for molecular properties."""

    cano_smiles: str
    num_heavy_atoms: int
    num_rings: int
    pfp: list[int]  # Pattern Fingerprint
