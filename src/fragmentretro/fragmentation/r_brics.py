# $Id$
#
#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#  Copyright (c) 2022, IBM LLC, Leili Zhang, Vasu Rao, Wendy Cornell
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#     * Neither the name of International Business Machine
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum, Nov 2008
# Updated by Leili Zhang, Vasu Rao, Jul 2022

"""Implementation of the BRICS algorithm from Degen et al. ChemMedChem *3* 1503-7 (2008)
*** Unpublished manuscript by Zhang et al. *** (to be updated)

"""

import copy
import random
import re
from collections.abc import Generator

from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions

from fragmentretro.exceptions import SmartsParsingError

type EnvironmentLabel = str
type SmartsPattern = str
type BondTypeSpecifier = str

type EnvironmentSmartsMap = dict[EnvironmentLabel, SmartsPattern]

type BondCleavageRule = tuple[EnvironmentLabel, EnvironmentLabel, BondTypeSpecifier]
type BondCleavageRuleGroup = list[BondCleavageRule]
type BondCleavageDefinitions = tuple[BondCleavageRuleGroup, ...]


type RDKitMol = Chem.Mol  # More explicit alias for rdkit.Chem.Mol
type RDKitReaction = Reactions.ChemicalReaction  # Alias for rdkit reaction object

type EnvironmentMolMatchersMap = dict[EnvironmentLabel, RDKitMol]

type BondMatcherRule = tuple[EnvironmentLabel, EnvironmentLabel, BondTypeSpecifier, RDKitMol]
type BondMatcherList = list[list[BondMatcherRule]]

# Type for a list of RDKit reaction objects
type ReactionList = list[RDKitReaction]
# Type for the grouped forward reactions
type ReactionGroups = tuple[ReactionList, ...]


# These are the definitions that will be applied to fragment molecules:
environs: EnvironmentSmartsMap = {
    "L1": "[C;D3]([#0,#6,#7,#8])(=O)",
    #
    # After some discussion, the L2 definitions ("N.pl3" in the original
    # paper) have been removed and incorporated into a (almost) general
    # purpose amine definition in L5 ("N.sp3" in the paper).
    #
    # The problem is one of consistency.
    #    Based on the original definitions you should get the following
    #    fragmentations:
    #      C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
    #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
    #    This difference just didn't make sense to us. By switching to
    #    the unified definition we end up with:
    #      C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
    #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
    #
    #'L2':'[N;!R;!D1;!$(N=*)]-;!@[#0,#6]',
    # this one turned out to be too tricky to define above, so we set it off
    # in its own definition:
    #'L2a':'[N;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]',
    "L3": "[O;D2]-;!@[#0,#6,#1]",
    "L4": "[C;!D1;!$(C=*)]-;!@[#6]",
    #'L5':'[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
    "L5": "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
    "L51": "[N;!R;!D1;$(N(!@[N,O]))]",  # Leili N-N, N-O, only breaks with carbon
    "L6": "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]",
    "L7a": "[C;D2,D3]-[#6]",
    "L7b": "[C;D2,D3]-[#6]",
    "#L8": "[C;!R;!D1]-;!@[#6]",  # Original L8
    "##L8": "[C;!R;!D1;!$(C!-*)]",
    "L8": "[C;!R;!D1;!$(C!-*);!$(C([H])([H])([H]))]",  # Leili fixed a problem when Hs are explicit
    "L81": "[C;!R;!D1;$(C(-[C,N,O,S])(=[N,S]))]",  # Leili amidine group
    #'L82':'[C;!R;!D1;$(C(-NH2)(=[O]))]', #Leili terminal COCH2, only for long chains
    "L9": "[RN,n;+0;$([RN,n](@[RC,RN,RO,RS,c,n,o,s])@[RC,RN,RO,RS,c,n,o,s])]",  # Leili generalized ring
    "#L9": "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]",  # Original L9
    "L10": "[N;R;$(N(@C(=O))@[C,N,O,S])]",
    "L11": "[S;D2](-;!@[#0,#6])",
    "L12": "[S;D4]([#6,#0])(=O)(=O)",
    "L12b": "[S;D4;!R](!@O)(!@O)",  # Leili, PDB-SMI conversion issue
    "L13": "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]",
    "L14": "[c;$(c(:[c,n,o,s]):[n,o,s])]",
    "L14b": "[RC;$([RC](@[RC,RN,RO,RS])@[RN,RO,RS])]",  # Leili generalized ring
    "L15": "[C;$(C(-;@C)-;@C)]",
    "L16": "[c;$(c(:c):c)]",
    "L16b": "[RC;$([RC](@[RC])@[RC])]",  # Leili generalized ring
    # Our Proposed Environments, by Zhang and Rao
    "L17": "[C](-C)(-C)(-C)",  # Aliphatic 1 - isobutyl/isopropyl
    "L18": "[R!#1;x3]",  # Ring 1
    "L19": "[R!#1;x2]",  # Ring 2
    "L182": "[R!#1;x3]",  # Ring 1b double bond - for assembly
    "L192": "[R!#1;x2]",  # Ring 2b double bond - for assembly
    "L20": "[CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",  # Aliphatic 2 - butyl, v2
    "L21": "[CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",  # Aliphatic 3 - butyl-CH3 or butyl-ring C or butyl-C-X, v2
    "L22": "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",  # Aliphatic 4 - octanyl, v2
    "L23": "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",  # Aliphatic 5 - octanyl-CH3 or octanyl-ring C or octonyl-C-X, v2
    "L30": "[C;D2]([#0,#6,#7,#8,#16])(#[N,C])",  # C#N C#C
    #'L51':'[N;!R;!D1;$(N(!@[N,O]))]', #duplication for record
    #'L81':'[C;!R;!D1;$(C(-[C,N,O,S])(=[N,S]))]', #duplication for record
}
reaction_defs: BondCleavageDefinitions = (
    # L1
    [
        ("1", "3", "-"),
        ("1", "5", "-"),
        ("1", "10", "-"),
    ],
    # L30
    [
        ("30", "30", "-"),
        ("30", "4", "-"),
        ("30", "5", "-"),
        ("30", "51", "-"),
        ("30", "6", "-"),
        ("30", "81", "-"),
        ("30", "9", "-"),
        ("30", "10", "-"),
        ("30", "11", "-"),
        ("30", "12", "-"),
        ("30", "12b", "-"),
        ("30", "13", "-"),
        ("30", "14", "-"),
        ("30", "14b", "-"),
        ("30", "15", "-"),
        ("30", "16", "-"),
        ("30", "16b", "-"),
    ],
    # L3
    [
        ("3", "4", "-"),
        ("3", "13", "-"),
        ("3", "14", "-"),
        ("3", "15", "-"),
        ("3", "16", "-"),
    ],
    # L4
    [
        ("4", "5", "-"),
        ("4", "11", "-"),
    ],
    # L5
    [
        ("5", "12", "-"),
        ("5", "14", "-"),
        ("5", "16", "-"),
        ("5", "13", "-"),
        ("5", "15", "-"),
    ],
    # L51 Leili
    [
        ("51", "1", "-"),
        ("51", "4", "-"),
        ("51", "12", "-"),
        ("51", "12b", "-"),
        ("51", "14", "-"),
        ("51", "16", "-"),
        ("51", "13", "-"),
        ("51", "15", "-"),
    ],
    # L6
    [
        ("6", "13", "-"),
        ("6", "14", "-"),
        ("6", "15", "-"),
        ("6", "16", "-"),
    ],
    # L7
    [
        ("7a", "7b", "="),
    ],
    # L8
    [
        ("8", "9", "-"),
        ("8", "10", "-"),
        ("8", "13", "-"),
        ("8", "14", "-"),
        ("8", "15", "-"),
        ("8", "16", "-"),
    ],
    # L81 Leili
    [
        ("81", "8", "-"),
        ("81", "9", "-"),
        ("81", "10", "-"),
        ("81", "13", "-"),
        ("81", "14", "-"),
        ("81", "15", "-"),
        ("81", "16", "-"),
    ],
    # L9
    [
        ("9", "13", "-"),  # not in Degen paper
        ("9", "14", "-"),  # not in Degen paper
        ("9", "15", "-"),
        ("9", "16", "-"),
    ],
    # L10
    [
        ("10", "13", "-"),
        ("10", "14", "-"),
        ("10", "15", "-"),
        ("10", "16", "-"),
    ],
    # L11
    [
        ("11", "13", "-"),
        ("11", "14", "-"),
        ("11", "15", "-"),
        ("11", "16", "-"),
    ],
    # L12
    # none left
    # L12b Leili
    [
        ("12b", "12b", "-"),
        ("12b", "5", "-"),
        ("12b", "4", "-"),
        ("12b", "13", "-"),
        ("12b", "14", "-"),
        ("12b", "15", "-"),
        ("12b", "16", "-"),
    ],
    # L13
    [
        ("13", "14", "-;@,!@"),
        ("13", "15", "-;!@"),
        ("13", "16", "-;@,!@"),
    ],
    # L14
    [
        ("14", "14", "-;@,!@"),
        ("14", "15", "-;@,!@"),
        ("14", "16", "-;@,!@"),
    ],
    # L14b
    [
        ("3", "14b", "-"),
        ("5", "14b", "-"),
        ("51", "14b", "-"),
        ("6", "14b", "-"),
        ("8", "14b", "-"),
        ("81", "14b", "-"),
        ("9", "14b", "-"),
        ("10", "14b", "-"),
        ("11", "14b", "-"),
        ("12b", "14b", "-"),
        ("13", "14b", "-"),
        ("14", "14b", "-"),
        ("14b", "14b", "-"),
        ("14b", "16", "-"),
        ("14b", "16b", "-"),
        ("14b", "15", "-"),
        ("14b", "17", "-"),
    ],
    # L15
    [
        ("15", "16", "-;@,!@"),
    ],
    # L16
    [
        ("16", "16", "-;@,!@"),  # not in Degen paper
    ],
    # L16b
    [
        ("3", "16b", "-"),
        ("5", "16b", "-"),
        ("51", "16b", "-"),
        ("6", "16b", "-"),
        ("8", "16b", "-"),
        ("81", "16b", "-"),
        ("9", "16b", "-"),
        ("10", "16b", "-"),
        ("11", "16b", "-"),
        ("12b", "16b", "-"),
        ("13", "16b", "-"),
        ("15", "16b", "-"),
        ("16", "16b", "-"),
        ("16b", "16b", "-"),
        ("17", "16b", "-"),
    ],
    # L17 Vasu and Leili (17-21)
    [
        ("17", "17", "-"),
        ("17", "16", "-"),
        ("17", "15", "-"),
        ("17", "14", "-"),
        ("17", "13", "-"),
        ("17", "12b", "-"),
        ("17", "12", "-"),
        ("17", "11", "-"),
        ("17", "10", "-"),
        ("17", "9", "-"),
        ("17", "8", "-"),
        ("17", "81", "-"),
        ("17", "51", "-"),
        ("17", "5", "-"),
    ],
    # L18
    [
        ("18", "19", "-;@"),
        ("182", "192", "=;@"),
    ],
    # L20
    [("20", "21", "-")],
    # L22
    [("22", "23", "-")],
)

reaction_defs_r: BondCleavageDefinitions = ([("20", "21", "-")],)


def initialize_env_matchers(environments: EnvironmentSmartsMap) -> EnvironmentMolMatchersMap:
    env_matchers = {}
    for env, smarts in environments.items():
        env_matchers[env] = Chem.MolFromSmarts(smarts)
    return env_matchers


def initialize_bond_matchers(
    environments: EnvironmentSmartsMap, reaction_definitions: BondCleavageDefinitions
) -> BondMatcherList:
    bond_matchers = []
    for rule_group in reaction_definitions:
        tmp = []
        for env_idx1, env_idx2, bond in rule_group:
            smarts1 = environments[f"L{env_idx1}"]
            smarts2 = environments[f"L{env_idx2}"]
            if "@" in bond:  # Leili
                patt = f"[$({smarts1})]{bond}[$({smarts2})]"
            else:
                # patt = "[$(%s)]%s;!@[$(%s)]" % (smarts1, bond, smarts2)
                patt = f"[$({smarts1})]{bond};!@[$({smarts2})]"
            # patt = '[$(%s)]%s;!@[$(%s)]'%(e1,bType,e2) #original
            patt = Chem.MolFromSmarts(patt)
            tmp.append((env_idx1, env_idx2, bond, patt))
        bond_matchers.append(tmp)
    return bond_matchers


def init_reactions(
    environments: EnvironmentSmartsMap, reaction_definitions: BondCleavageDefinitions
) -> tuple[list[list[SmartsPattern]], ReactionGroups, ReactionList]:
    rule_groups = copy.deepcopy(reaction_definitions)
    smarts_groups = []
    for rule_group in rule_groups:
        new_rule_group = []
        for rule in rule_group:
            env_idx1, env_idx2, bond = rule
            smarts1 = environs["L" + env_idx1]
            smarts2 = environs["L" + env_idx2]
            g1 = re.sub("[a-z,A-Z]", "", env_idx1)
            g2 = re.sub("[a-z,A-Z]", "", env_idx2)
            if "@" not in bond:  # Leili
                # if 1 == 1:
                # sma = "[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]" % (smarts1, bond, smarts2, g1, g2)
                sma = f"[$({smarts1})]{bond};!@[$({smarts2})]>>[{g1}*]-[*:1].[{g2}*]-[*:2]"
            else:
                # sma = "[$(%s):1]%s[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]" % (smarts1, bond, smarts2, g1, g2)
                sma = f"[$({smarts1})]{bond}[$({smarts2})]>>[{g1}*]-[*:1].[{g2}*]-[*:2]"
            new_rule_group.append(sma)
            # sma='[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]'%(r1,bnd,r2,g1,g2) #original
            # gp[j] =sma
        smarts_groups.append(new_rule_group)

    for smarts_group in smarts_groups:
        for smarts in smarts_group:
            try:
                t = Reactions.ReactionFromSmarts(smarts)
                t.Initialize()
            except Exception as e:
                raise SmartsParsingError(f"Failed to parse SMARTS pattern: {smarts}") from e

    reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smarts_groups])
    reverse_reactions = []
    for smarts_group in smarts_groups:
        for smarts in smarts_group:
            reactants, products = smarts.split(">>")
            smarts = f"{products}>>{reactants}"
            rxn = Reactions.ReactionFromSmarts(smarts)
            labels = re.findall(r"\[([0-9]+?)\*\]", products)
            rxn._matchers = [Chem.MolFromSmiles(f"[{x}*]") for x in labels]
            reverse_reactions.append(rxn)

    return smarts_groups, reactions, reverse_reactions


smarts_groups, reactions, REVERSE_REACTIONS = init_reactions(environs, reaction_defs)
smarts_groups_r, reactions_r, REVERSE_REACTIONS_R = init_reactions(environs, reaction_defs_r)

ENV_MATCHERS = initialize_env_matchers(environs)
# AM: rBRICS had env_matchers_r, but because env matcher does not depend on reaction definitions,
# that env_matchers_r was identical to env_matchers

BOND_MATCHERS = initialize_bond_matchers(environs, reaction_defs)
BOND_MATCHERS_R = initialize_bond_matchers(environs, reaction_defs_r)


def find_brics_bonds(
    mol: Chem.Mol, randomizeOrder: bool = False, silent: bool = True
) -> Generator[tuple[tuple[int, int], tuple[EnvironmentLabel, EnvironmentLabel]], None, None]:
    """returns the bonds in a molecule that BRICS would cleave

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(find_brics_bonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4'))]

    a more complicated case:
    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> res = list(find_brics_bonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    we can also randomize the order of the results:
    >>> random.seed(23)
    >>> res = list(find_brics_bonds(m,randomizeOrder=True))
    >>> sorted(res)
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    Note that this is a generator function :
    >>> res = find_brics_bonds(m)
    >>> res
    <generator object ...>
    >>> res.next()
    ((3, 2), ('3', '4'))

    >>> m = Chem.MolFromSmiles('CC=CC')
    >>> res = list(find_brics_bonds(m))
    >>> sorted(res)
    [((1, 2), ('7', '7'))]

    make sure we don't match ring bonds:
    >>> m = Chem.MolFromSmiles('O=C1NCCC1')
    >>> list(find_brics_bonds(m))
    []

    another nice one, make sure environment 8 doesn't match something connected
    to a ring atom:
    >>> m = Chem.MolFromSmiles('CC1(C)CCCCC1')
    >>> list(FindrBRICSBonds(m))
    []

    """
    letter = re.compile("[a-z,A-Z]")
    indices = list(range(len(BOND_MATCHERS)))
    bonds_done = set()
    if randomizeOrder:
        random.shuffle(indices)

    env_matches = {}
    for env, substruct_mol in ENV_MATCHERS.items():  # Leili
        env_matches[env] = mol.HasSubstructMatch(substruct_mol)
    for group_idx in indices:
        if randomizeOrder:
            matchers_list = BOND_MATCHERS[group_idx][:]
            random.shuffle(matchers_list)
        else:
            matchers_list = BOND_MATCHERS[group_idx]
        for i1, i2, _, substruct_mol in matchers_list:
            if not env_matches["L" + i1] or not env_matches["L" + i2]:
                continue
            matches = mol.GetSubstructMatches(substruct_mol)
            i1 = letter.sub("", i1)
            i2 = letter.sub("", i2)
            for match in matches:
                if match not in bonds_done and (match[1], match[0]) not in bonds_done:
                    bonds_done.add(match)
                    yield (((match[0], match[1]), (i1, i2)))


# Leili
def find_r_brics_bonds(
    mol: Chem.Mol, randomizeOrder: bool = False, silent: bool = True
) -> Generator[tuple[tuple[int, int], tuple[EnvironmentLabel, EnvironmentLabel]], None, None]:
    letter = re.compile("[a-z,A-Z]")
    indices = list(range(len(BOND_MATCHERS_R)))

    bonds_done = set()
    if randomizeOrder:
        random.shuffle(indices)

    env_matches = {}
    for env, substruct_mol in ENV_MATCHERS.items():  # Leili
        env_matches[env] = mol.HasSubstructMatch(substruct_mol)
    for group_idx in indices:
        if randomizeOrder:
            matchers_list = BOND_MATCHERS_R[group_idx][:]
            random.shuffle(matchers_list)
        else:
            matchers_list = BOND_MATCHERS_R[group_idx]
        for i1, i2, _, substruct_mol in matchers_list:
            if not env_matches["L" + i1] or not env_matches["L" + i2]:
                continue
            matches = mol.GetSubstructMatches(substruct_mol)
            i1 = letter.sub("", i1)
            i2 = letter.sub("", i2)
            for match in matches:
                if match not in bonds_done and (match[1], match[0]) not in bonds_done:
                    bonds_done.add(match)
                    yield (((match[0], match[1]), (i1, i2)))


def break_r_brics_bonds(
    mol: Chem.Mol,
    bonds: list[tuple[tuple[int, int], tuple[EnvironmentLabel, EnvironmentLabel]]] | None = None,
    sanitize: bool = True,
    silent: bool = True,
) -> Chem.Mol:
    """breaks the BRICS bonds in a molecule and returns the results

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=break_r_brics_bonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[3*]O[3*].[4*]CC.[4*]CCC'

    a more complicated case:
    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> m2=break_r_brics_bonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O.[16*]c1ccccc1'


    can also specify a limited set of bonds to work with:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2 = break_r_brics_bonds(m,[((3, 2), ('3', '4'))])
    >>> Chem.MolToSmiles(m2,True)
    '[3*]OCC.[4*]CCC'

    this can be used as an alternate approach for doing a BRICS decomposition by
    following break_r_brics_bonds with a call to Chem.GetMolFrags:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=break_r_brics_bonds(m)
    >>> frags = Chem.GetMolFrags(m2,asMols=True)
    >>> [Chem.MolToSmiles(x,True) for x in frags]
    ['[4*]CCC', '[3*]O[3*]', '[4*]CC']

    """
    if not bonds:
        # bonds = FindrBRICSBonds(mol)
        res = Chem.FragmentOnBRICSBonds(mol)
        if sanitize:
            Chem.SanitizeMol(res)
        return res
    eMol = Chem.EditableMol(mol)
    nAts = mol.GetNumAtoms()

    dummyPositions = []
    for indices, dummyTypes in bonds:
        ia, ib = indices
        obond = mol.GetBondBetweenAtoms(ia, ib)
        bondType = obond.GetBondType()
        eMol.RemoveBond(ia, ib)

        da, db = dummyTypes
        atoma = Chem.Atom(0)
        atoma.SetIsotope(int(da))
        atoma.SetNoImplicit(True)
        idxa = nAts
        nAts += 1
        eMol.AddAtom(atoma)
        eMol.AddBond(ia, idxa, bondType)

        atomb = Chem.Atom(0)
        atomb.SetIsotope(int(db))
        atomb.SetNoImplicit(True)
        idxb = nAts
        nAts += 1
        eMol.AddAtom(atomb)
        eMol.AddBond(ib, idxb, bondType)
        if mol.GetNumConformers():
            dummyPositions.append((idxa, ib))
            dummyPositions.append((idxb, ia))
    res = eMol.GetMol()
    if sanitize:
        Chem.SanitizeMol(res)
    if mol.GetNumConformers():
        for conf in mol.GetConformers():
            resConf = res.GetConformer(conf.GetId())
            for ia, pa in dummyPositions:
                resConf.SetAtomPosition(ia, conf.GetAtomPosition(pa))
    return res


def r_brics_decompose(
    mol: Chem.Mol,
    allNodes: set[str] | None = None,
    minFragmentSize: int = 1,
    onlyUseReactions: list[tuple[int, int]] | None = None,
    silent: bool = True,
    keepNonLeafNodes: bool = False,
    singlePass: bool = False,
    returnMols: bool = False,
) -> list[Chem.Mol]:
    """returns the BRICS decomposition for a molecule

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(r_brics_decompose(m))
    >>> sorted(res)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    >>> res = r_brics_decompose(m,returnMols=True)
    >>> res[0]
    <rdkit.Chem.rdchem.Mol object ...>
    >>> smis = [Chem.MolToSmiles(x,True) for x in res]
    >>> sorted(smis)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    nexavar, an example from the paper (corrected):
    >>> m = Chem.MolFromSmiles('CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
    >>> res = list(r_brics_decompose(m))
    >>> sorted(res)
    ['[1*]C([1*])=O', '[1*]C([6*])=O', '[14*]c1cc([16*])ccn1', '[16*]c1ccc(Cl)c([16*])c1', '[16*]c1ccc([16*])cc1', '[3*]O[3*]', '[5*]NC', '[5*]N[5*]', '[8*]C(F)(F)F']

    it's also possible to keep pieces that haven't been fully decomposed:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(r_brics_decompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[3*]O[3*]', '[4*]CC', '[4*]CCC']

    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(r_brics_decompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[16*]c1cccc([16*])c1', '[3*]OCCC', '[3*]OC[8*]', '[3*]OCc1cccc(-c2ccccn2)c1', '[3*]OCc1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]', '[4*]Cc1cccc(-c2ccccn2)c1', '[4*]Cc1cccc([16*])c1', '[8*]COCCC']

    or to only do a single pass of decomposition:
    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(r_brics_decompose(m,singlePass=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[3*]OCCC', '[3*]OCc1cccc(-c2ccccn2)c1', '[4*]CCC', '[4*]Cc1cccc(-c2ccccn2)c1', '[8*]COCCC']

    setting a minimum size for the fragments:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(r_brics_decompose(m,keepNonLeafNodes=True,minFragmentSize=2))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(r_brics_decompose(m,keepNonLeafNodes=True,minFragmentSize=3))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[4*]CCC']
    >>> res = list(r_brics_decompose(m,minFragmentSize=2))
    >>> sorted(res)
    ['[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']


    """
    global reactions
    mSmi = Chem.MolToSmiles(mol, 1)

    if allNodes is None:
        allNodes = set()

    if mSmi in allNodes:
        return []

    activePool: dict[str, Chem.Mol] = {mSmi: mol}
    allNodes.add(mSmi)
    foundMols: dict[str, Chem.Mol] = {mSmi: mol}
    for gpIdx, reactionGp in enumerate(reactions):
        newPool: dict[str, Chem.Mol] = {}
        while activePool:
            matched = False
            nSmi = list(activePool.keys())[0]
            mol = activePool.pop(nSmi)
            for rxnIdx, reaction in enumerate(reactionGp):
                if onlyUseReactions and (gpIdx, rxnIdx) not in onlyUseReactions:
                    continue
                # if not silent:
                #     print("--------")
                #     print(smartsGps[gpIdx][rxnIdx])
                ps = reaction.RunReactants((mol,))
                if ps:
                    if not silent:
                        print(nSmi, "->", len(ps), "products")
                    for prodSeq in ps:
                        seqOk = True
                        # we want to disqualify small fragments, so sort the product sequence by size
                        prodSeq = [(prod.GetNumAtoms(onlyExplicit=True), prod) for prod in prodSeq]
                        print(prodSeq)
                        prodSeq.sort(key=lambda y: y[0])
                        for nats, prod in prodSeq:
                            try:
                                Chem.SanitizeMol(prod)
                            except Exception:
                                continue
                            pSmi = Chem.MolToSmiles(prod, 1)
                            if minFragmentSize > 0:
                                nDummies = pSmi.count("*")
                                if nats - nDummies < minFragmentSize:
                                    seqOk = False
                                    break
                            prod.pSmi = pSmi
                        if seqOk:
                            matched = True
                            for _, prod in prodSeq:
                                pSmi = prod.pSmi
                                # print '\t',nats,pSmi
                                if pSmi not in allNodes:
                                    if not singlePass:
                                        activePool[pSmi] = prod
                                    allNodes.add(pSmi)
                                    foundMols[pSmi] = prod
            if singlePass or keepNonLeafNodes or not matched:
                newPool[nSmi] = mol
        activePool = newPool
    if not (singlePass or keepNonLeafNodes):
        if not returnMols:
            return list(activePool.keys())
        else:
            return list(activePool.values())
    else:
        if not returnMols:
            return list(allNodes)
        else:
            return list(foundMols.values())


dummyPattern = Chem.MolFromSmiles("[*]")


def BRICSBuild(
    fragments: list[Chem.Mol],
    onlyCompleteMols: bool = True,
    seeds: list[Chem.Mol] | None = None,
    uniquify: bool = True,
    scrambleReagents: bool = True,
    maxDepth: int = 3,
) -> Generator[Chem.Mol, None, None]:
    seen = set()
    if not seeds:
        seeds = list(fragments)
    if scrambleReagents:
        seeds = list(seeds)
        random.shuffle(seeds)
    if scrambleReagents:
        tempReactions = list(REVERSE_REACTIONS)
        random.shuffle(tempReactions)
    else:
        tempReactions = REVERSE_REACTIONS
    for seed in seeds:
        seedIsR1 = False
        seedIsR2 = False
        nextSteps = []
        for rxn in tempReactions:
            if seed.HasSubstructMatch(rxn._matchers[0]):
                seedIsR1 = True
            if seed.HasSubstructMatch(rxn._matchers[1]):
                seedIsR2 = True
            for fragment in fragments:
                ps = None
                if fragment.HasSubstructMatch(rxn._matchers[0]) and seedIsR2:
                    ps = rxn.RunReactants((fragment, seed))
                if fragment.HasSubstructMatch(rxn._matchers[1]) and seedIsR1:
                    ps = rxn.RunReactants((seed, fragment))
                if ps:
                    for p in ps:
                        if uniquify:
                            pSmi = Chem.MolToSmiles(p[0], True)
                            if pSmi in seen:
                                continue
                            else:
                                seen.add(pSmi)
                        if p[0].HasSubstructMatch(dummyPattern):
                            nextSteps.append(p[0])
                            if not onlyCompleteMols:
                                yield p[0]
                        else:
                            yield p[0]
        if nextSteps and maxDepth > 0:
            for p in BRICSBuild(
                fragments, onlyCompleteMols=onlyCompleteMols, seeds=nextSteps, uniquify=uniquify, maxDepth=maxDepth - 1
            ):
                if uniquify:
                    pSmi = Chem.MolToSmiles(p, True)
                    if pSmi in seen:
                        continue
                    else:
                        seen.add(pSmi)
                yield p


#
# Leili
# Iteratively breaks all aliphatic chains using linkage L20-L21, just in case chain is too long and L20-L23 aren't enough
# Not the most efficient code yet
def reBRICS(fragments: list[Chem.Mol]) -> list[Chem.Mol]:
    oldfragments = fragments
    breakable = [1] * len(oldfragments)
    breakout = sum(breakable)
    iteri = 0
    while breakout > 0:
        newfrags: list[Chem.Mol] = []
        newbreakable: list[int] = []
        for frag in oldfragments:
            heavy = frag.GetNumHeavyAtoms()
            if heavy > 5:
                if len(frag.GetSubstructMatch(Chem.MolFromSmiles("CCCCCC"))) > 0:
                    secondbonds = list(find_r_brics_bonds(frag))  # only breaks aliphatic chains
                    secondpieces = break_r_brics_bonds(frag, secondbonds)
                    secondfrags = Chem.GetMolFrags(secondpieces, asMols=True)
                    if len(secondfrags) > 1:
                        newfrags = newfrags + secondfrags
                        newbreakable = newbreakable + [1] * len(secondfrags)
                    else:
                        newfrags = newfrags + secondfrags
                        newbreakable.append(0)  # 1 frag->1 frag, no longer breakable, criteria #1
                else:
                    newfrags.append(frag)
                    newbreakable.append(0)  # frag doesn't contain CCCCCC chain, no longer breakable, criteria #2
            else:
                newfrags.append(frag)
                newbreakable.append(0)  # frag contains less than 5 heavy atoms, no longer breakable, criteria #3
        oldfragments = newfrags
        breakable = newbreakable
        breakout = sum(newbreakable)
        iteri += 1
        if iteri > 100:  # more than 100 iterations, criteria #4
            print("Too many iterations. Not trying to break a polymer, are we? :)")
            breakout = 0
    return newfrags
