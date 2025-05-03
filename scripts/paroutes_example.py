import pickle
from pathlib import Path

from tqdm import tqdm

from fragmentretro.fragmenter import BRICSFragmenter  # or use rBRICSFragmenter
from fragmentretro.retrosynthesis import Retrosynthesis
from fragmentretro.solutions import RetrosynthesisSolution
from fragmentretro.utils.logging_config import logger

DATA_PATH = Path(__name__).parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
EVAL_PATH = DATA_PATH / "evaluations"
PRECOMPUTE_PATH = DATA_PATH / "precompute"


####### Change this
with open(PAROUTES_PATH / "n1-stock.txt") as f:
    n1_stock = [line.strip() for line in f.readlines()]

with open(PAROUTES_PATH / "n1-targets.txt", "rb") as f:
    n1_targets = pickle.load(f)

# with open(PAROUTES_PATH / "n5-stock.txt", "r") as f:
#     n5_stock = [line.strip() for line in f.readlines()]

# with open(PAROUTES_PATH / "n5-targets.txt", "rb") as f:
#     n5_targets = pickle.load(f)


targets = n1_targets
stock = n1_stock
SAVED_PATH = EVAL_PATH / "save_path.pkl"

original_BBs = None
JSON_PRECOMPUTE_PATH = PRECOMPUTE_PATH / "n1_stock_properties.json"

# original_BBs = set(stock)
# JSON_PRECOMPUTE_PATH = None
#######


all_retro_tool = []
all_solution_tool = []
solved_count = 0
for i in tqdm(range(len(targets))):
    target = targets[i]
    fragmenter = BRICSFragmenter(target)
    retro_tool = Retrosynthesis(fragmenter, original_BBs=original_BBs, mol_properties_path=JSON_PRECOMPUTE_PATH)
    retro_tool.fragment_retrosynthesis()  # run retrosynthesis
    retro_solution = RetrosynthesisSolution(retro_tool)
    retro_solution.fill_solutions()
    logger.info(f"Target {i + 1}, # of solutions: {len(retro_solution.solutions)}")
    all_retro_tool.append(retro_tool)
    all_solution_tool.append(retro_solution)
    if len(retro_solution.solutions) > 0:
        solved_count += 1
    logger.info(f"Current Solved Rate: {solved_count / (i + 1)}")

logger.info(f"Final Solved Rate: {solved_count / len(targets)}")

with open(SAVED_PATH, "wb") as f:
    pickle.dump([all_retro_tool, all_solution_tool], f)
