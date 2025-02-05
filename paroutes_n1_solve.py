import pickle
from pathlib import Path

from tqdm import tqdm

from FragmentRetro.fragmenter import BRICSFragmenter
from FragmentRetro.logging_config import logger
from FragmentRetro.retrosynthesis import Retrosynthesis
from FragmentRetro.solutions import RetrosynthesisSolution

DATA_PATH = Path(__name__).parent / "data"
PAROUTES_PATH = DATA_PATH / "paroutes"
EVAL_PATH = DATA_PATH / "evaluations"


####### Change this
with open(PAROUTES_PATH / "n1-stock.txt", "r") as f:
    n1_stock = [line.strip() for line in f.readlines()]

with open(PAROUTES_PATH / "n1_shuffled_seed42_n500.pkl", "rb") as f:
    n1_data = pickle.load(f)

n1_products, n1_routes = n1_data[0], n1_data[2]

targets = n1_products
stock = n1_stock
SAVED_PATH = EVAL_PATH / "eval_n1_500_seed42.pkl"
#######


all_retro_tool = []
all_solution_tool = []
solved_count = 0
for i in tqdm(range(len(targets))):
    target = targets[i]
    fragmenter = BRICSFragmenter(target)
    retro_tool = Retrosynthesis(fragmenter, set(stock))
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
