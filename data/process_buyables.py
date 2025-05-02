import gzip
import json
import shutil
from pathlib import Path

from tqdm import tqdm

from FragmentRetro.utils.filter_compound import precompute_properties
from FragmentRetro.utils.helpers import canonicalize_smiles
from FragmentRetro.utils.logging_config import logger

DATA_PATH = Path(__file__).parent
BUYABLES_PATH = DATA_PATH / "buyables"
PRECOMPUTE_PATH = DATA_PATH / "precompute"

if __name__ == "__main__":
    logger.info("Decompressing buyables data...")
    with gzip.open(BUYABLES_PATH / "buyables_all.json.gz", "rb") as f_in:
        with open(BUYABLES_PATH / "buyables_all.json", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    logger.info("Loading buyables data...")
    with open(BUYABLES_PATH / "buyables_all.json", "r") as f:
        buyables = json.load(f)
    buyables_smiles = [buyable["smiles"] for buyable in buyables]

    logger.info("Canonicalizing SMILES strings...")
    canonicalized_smiles = []
    for smiles in tqdm(buyables_smiles):
        try:
            canonicalized_smiles.append(canonicalize_smiles(smiles))
        except ValueError as e:
            logger.error(f"Error canonicalizing SMILES '{smiles}': {e} during canonicalizing buyables")

    logger.info("Saving unique canonicalized SMILES strings...")
    unique_smiles = list(set(canonicalized_smiles))
    with open(BUYABLES_PATH / "buyables-stock.txt", "w") as f:
        for smiles in unique_smiles:
            f.write(smiles + "\n")

    logger.info("Precomputing properties for buyables...")
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(unique_smiles, PRECOMPUTE_PATH / "buyables_stock_properties.json")
