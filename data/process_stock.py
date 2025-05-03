from pathlib import Path

from FragmentRetro.utils.filter_compound import precompute_properties

DATA_PATH = Path(__file__).parent
PAROUTES_PATH = DATA_PATH / "paroutes"
PRECOMPUTE_PATH = DATA_PATH / "precompute"


if __name__ == "__main__":
    with open(PAROUTES_PATH / "n1-stock.txt") as f:
        n1_stock = [line.strip() for line in f.readlines()]
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(n1_stock, PRECOMPUTE_PATH / "n1_stock_properties.json")
    # precompute_properties(n1_stock, PRECOMPUTE_PATH / "n1_stock_properties_fp1024.json", fpSize=1024)

    with open(PAROUTES_PATH / "n5-stock.txt") as f:
        n5_stock = [line.strip() for line in f.readlines()]
    PRECOMPUTE_PATH.mkdir(parents=True, exist_ok=True)
    precompute_properties(n5_stock, PRECOMPUTE_PATH / "n5_stock_properties.json")
    # precompute_properties(n5_stock, PRECOMPUTE_PATH / "n5_stock_properties_fp1024.json", fpSize=1024)
