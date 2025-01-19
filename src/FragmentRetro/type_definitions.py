"""Type definitions for FragmentRetro package."""

from typing import List, Tuple

from typing_extensions import TypeAlias

# Bond type (from BRICS): ((atom1_idx, atom2_idx), (type1, type2))
# Example: ((3, 2), ('3', '4'))
BondType: TypeAlias = Tuple[Tuple[int, int], Tuple[str, str]]
BondsType: TypeAlias = List[BondType]

# Atom mapping type: tuple of atom indices
# Example: (0, 1, 2, 7)
AtomMappingType: TypeAlias = Tuple[int, ...]
AtomMappingsType: TypeAlias = List[AtomMappingType]
