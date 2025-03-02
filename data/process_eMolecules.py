from pathlib import Path

import pandas as pd
from tqdm import tqdm

from FragmentRetro.utils.filter_compound import precompute_properties
from FragmentRetro.utils.helpers import canonicalize_smiles
from FragmentRetro.utils.logging_config import logger

DATA_PATH = Path(__file__).parent
EMOLECULES_PATH = DATA_PATH / "eMolecules"
PRECOMPUTE_PATH = DATA_PATH / "precompute"

if __name__ == "__main__":
    logger.info("Loading eMolecules csv...")
    emol_df = pd.read_csv(EMOLECULES_PATH / "origin_dict.csv", index_col=0)
    emol_smiles = emol_df["mol"].tolist()
    logger.info(f"Number of eMolecules smiles: {len(emol_smiles)}")
    del emol_df

    logger.info("Canonicalizing SMILES strings...")
    canonicalized_smiles = []
    for smiles in tqdm(emol_smiles):
        try:
            canonicalized_smiles.append(canonicalize_smiles(smiles))
        except ValueError as e:
            logger.error(f"Error canonicalizing SMILES '{smiles}': {e} during canonicalizing buyables")

    logger.info("Saving unique canonicalized SMILES strings...")
    unique_smiles = list(set(canonicalized_smiles))
    logger.info(f"Number of unique eMolecules smiles: {len(unique_smiles)}")
    with open(EMOLECULES_PATH / "eMolecules.txt", "w") as f:
        for smiles in unique_smiles:
            f.write(smiles + "\n")

    logger.info("Precomputing properties for buyables...")
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(unique_smiles, PRECOMPUTE_PATH / "eMolecules_stock_properties.json")
