import itertools
from abc import ABC, abstractmethod
from typing import cast

import matplotlib.pyplot as plt
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Mol

from FragmentRetro.utils.logging_config import logger
from FragmentRetro.utils.type_definitions import AtomMappingType, BondType, CombType


class Fragmenter(ABC):
    def __init__(self, smiles: str, clearAromaticFlags: bool = False) -> None:
        self.original_smiles: str = smiles
        self.original_mol: Mol = Chem.MolFromSmiles(smiles)
        if clearAromaticFlags:
            # can break aromatic bonds
            Chem.Kekulize(self.original_mol, clearAromaticFlags=True)
        self.fragmentation_bonds: list[BondType] = self._find_fragmentation_bonds(self.original_mol)
        self.broken_mol: Mol = self._break_bonds(self.original_mol, self.fragmentation_bonds)
        self.atom_mappings: list[AtomMappingType] = []
        self.fragment_graph: nx.Graph = self._build_fragment_graph()
        self.num_fragments: int = len(self.fragment_graph.nodes())

    @abstractmethod
    def _find_fragmentation_bonds(self, mol: Mol) -> list[BondType]:
        pass

    @abstractmethod
    def _break_bonds(self, mol: Mol, bonds: list[BondType]) -> Mol:
        pass

    def _build_fragment_graph(self) -> nx.Graph:
        """
        Build graph representing fragment connectivity.
        Nodes are fragments with their SMILES and atom mappings.
        Edges contain bond type information.
        Uses Graph to support only one edge between fragments.

        Returns:
            NetworkX Graph representing fragment connectivity
        """
        G = nx.Graph()

        mol_fragments: tuple[Mol] = Chem.GetMolFrags(
            self.broken_mol, asMols=True, fragsMolAtomMapping=self.atom_mappings
        )

        atom_to_frag = {}
        for frag_id, atom_indices in enumerate(self.atom_mappings):
            for atom_idx in atom_indices:
                atom_to_frag[atom_idx] = frag_id

        for i, fragment in enumerate(mol_fragments):
            G.add_node(i, smiles=Chem.MolToSmiles(fragment), atom_indices=self.atom_mappings[i])

        edge_index = 0
        for (atom1, atom2), (type1, type2) in self.fragmentation_bonds:
            frag1 = atom_to_frag.get(atom1)
            frag2 = atom_to_frag.get(atom2)

            if frag1 != frag2 and frag1 is not None and frag2 is not None:
                if not G.has_edge(frag1, frag2):
                    G.add_edge(
                        frag1,
                        frag2,
                        bond_type=(type1, type2),
                        atoms=(atom1, atom2),
                        edge_index=edge_index,
                    )
                    edge_index += 1

        return G

    def _get_initial_fragments(self) -> list[str]:
        """
        Retrieve the initial fragments as SMILES strings from the fragment graph.

        Returns:
            List of SMILES strings representing the initial fragments.
        """
        return [self.fragment_graph.nodes[node]["smiles"] for node in self.fragment_graph.nodes()]

    def visualize(
        self,
        figsize: tuple[float, float] = (10.0, 10.0),
        with_indices: bool = False,
    ) -> None:
        """
        Visualize the fragment graph and optionally print detailed information.
        Handles multiple edges between nodes.
        """
        pos = nx.spring_layout(self.fragment_graph)
        plt.figure(figsize=figsize)

        # Draw nodes
        nx.draw_networkx_nodes(self.fragment_graph, pos, node_color="lightblue", node_size=2000)

        # Draw edges without curves since we only have single edges
        for edge in self.fragment_graph.edges(data=True):
            u, v, data = edge
            nx.draw_networkx_edges(
                self.fragment_graph,
                pos,
                edgelist=[(u, v)],
            )

        # Add node labels
        labels = {
            node: f"Node {node}: {data['smiles']}" if with_indices else data["smiles"]
            for node, data in self.fragment_graph.nodes(data=True)
        }
        nx.draw_networkx_labels(self.fragment_graph, pos, labels)

        # Add edge labels
        edge_labels = {
            (u, v): f"Edge {data['edge_index']} ({u}-{v}): {data['bond_type']}"
            if with_indices
            else f"{data['bond_type']}"
            for u, v, data in self.fragment_graph.edges(data=True)
        }

        for (u, v), label in edge_labels.items():
            x = (pos[u][0] + pos[v][0]) / 2
            y = (pos[u][1] + pos[v][1]) / 2
            plt.text(
                x,
                y,
                label,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
                horizontalalignment="center",
                verticalalignment="center",
            )

        plt.title("Fragment Connectivity Graph")
        plt.axis("off")
        plt.show()

        logger.info("[Fragmenter] \nNode data:")
        for node in self.fragment_graph.nodes():
            logger.info(f"[Fragmenter] \nNode {node}:")
            logger.info(f"[Fragmenter] SMILES: {self.fragment_graph.nodes[node]['smiles']}")
            logger.info(f"[Fragmenter] Atom indices: {self.fragment_graph.nodes[node]['atom_indices']}")

        logger.info("[Fragmenter] \nEdge data:")
        for u, v, data in self.fragment_graph.edges(data=True):
            logger.info(f"[Fragmenter] \nEdge {data['edge_index']} ({u}-{v}):")
            logger.info(f"[Fragmenter] Bond type: {data['bond_type']}")
            logger.info(f"[Fragmenter] Atoms: {data['atoms']}")

    def get_length_n_combinations(self, n: int) -> set[CombType]:
        """
        Get all unique combinations of n fragments in the fragment graph.

        Args:
            n: Length of combinations to find

        Returns:
            Set of unique combinations as sorted lists of node IDs
        """
        all_combinations = set()

        def dfs(path: list[int], candidates: set[int]) -> None:
            if len(path) == n:
                sorted_path = cast(CombType, tuple(sorted(path)))
                all_combinations.add(sorted_path)
                return

            for node in candidates:
                if node not in path:
                    new_candidates = candidates | set(self.fragment_graph.neighbors(node)) - set(path)
                    dfs(path + [node], new_candidates)

        for node in self.fragment_graph.nodes:
            dfs([node], set(self.fragment_graph.neighbors(node)))

        return all_combinations

    def check_connected_subgraph(self, combination: CombType) -> bool:
        """Check if the combination is a connected subgraph of the original molecule.

        Args:
            check_subgraph: Whether to check if the combination is a connected subgraph of the original molecule.

        Returns:
            True if the combination is a connected subgraph, False otherwise.

        """
        # check if the combination is a connected subgraph
        subgraph = self.fragment_graph.subgraph(combination)
        if not nx.is_connected(subgraph):
            return False
        return True

    def get_combination_smiles(self, combination: CombType) -> str:
        """
        Get the fragment smiles given one combination.

        Args:
            combination: A combination as a sorted list of node IDs.

        Returns:
            A SMILES string representing the fragment combination.
        """
        if len(combination) == 0:
            return ""
        elif len(combination) == 1:
            return cast(str, self.fragment_graph.nodes[combination[0]]["smiles"])
        elif len(combination) == self.num_fragments:
            return self.original_smiles

        # remove the bonds that are within the fragment combination
        bonds_to_break: list[BondType] = self.fragmentation_bonds.copy()
        for pair in itertools.combinations(combination, 2):
            edge_data = self.fragment_graph.get_edge_data(*pair)
            if edge_data:
                bonds_to_break.remove((edge_data["atoms"], edge_data["bond_type"]))
        # break the bonds
        comb_broken_mol = self._break_bonds(self.original_mol, bonds_to_break)
        comb_atom_mappings: list[AtomMappingType] = []
        mol_fragments: tuple[Mol] = Chem.GetMolFrags(
            comb_broken_mol, asMols=True, fragsMolAtomMapping=comb_atom_mappings
        )
        # get the smiles with the atom mapping that contains the first atom of the combination
        # first atom is chosen since it is always not the bond break (i.e. any or *) atom
        for i, mol in enumerate(mol_fragments):
            if self.atom_mappings[combination[0]][0] in comb_atom_mappings[i]:
                comb_smiles = Chem.MolToSmiles(mol)
                break
        return cast(str, comb_smiles)
