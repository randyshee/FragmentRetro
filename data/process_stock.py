from FragmentRetro.utils.filter_compounds import precompute_properties

from pathlib import Path

DATA_PATH = Path(__name__).parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
PRECOMPUTE_PATH = DATA_PATH / "precompute"


if __name__ == '__main__':
    with open(PAROUTES_PATH / "n1-stock.txt", 'r') as f:
        n1_stock = [line.strip() for line in f.readlines()]
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(n1_stock, PRECOMPUTE_PATH / "n1_stock_properties.json")

    with open(PAROUTES_PATH / "n5-stock.txt", 'r') as f:
        n5_stock = [line.strip() for line in f.readlines()]
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(n5_stock, PRECOMPUTE_PATH / "n5_stock_properties.json")
