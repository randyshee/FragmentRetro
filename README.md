# FragmentRetro

## Development

If you need other packages to run your features, you must specify the dependency in [pyproject.toml](./pyproject.toml) and the version used.

### Set Up

```bash
uv venv --python 3.11.4
```

which will create a new virtual environment (in the folder `venv` which is already gitignored).

Activate the environment:

```bash
source .venv/bin/activate
```

and install all dependencies (specified in [pyproject.toml](./pyproject.toml)):

```bash
uv pip install -e ".[dev]"
```

If you need to use jupyter notebook:

```bash
uv pip install ipykernel
```

### Cloning & Setting Pre-commits

1. Clone the repository.
2. Inside of the repository, run `pre-commit install`
3. Create a new branch `git branch nameofthebranch` and switch to it `git checkout nameofthebranch`.
4. Whenever you want to push the commits to the repo, do `git push -u origin nameofthebranch`

#### Pre-Commits

After `pre-commit install`, each time you try to make a commit, the system will automatically check for:

- linting errors (using [ruff](https://github.com/astral-sh/ruff))
- formatting (using [black](https://github.com/psf/black))
- import sorting (using [isort](https://github.com/pycqa/isort/)) - this is arguably the least necessary thing, but once you have it as a pre-commit hook, you forget it even exists.

> mypy from pre-commit hooks currently disabled just so that you can push commits even with failing mypy so that we can discuss it on GitHub and see how it can be fixed

If any of the checks above (ruff/black/isort) fail, the system will try to fix them automatically, so in most cases you'll only need to run `git add .` and try to commit again. If you want run the checks manually, you can do:

- `ruff check`. If there are errors, most of the time they can be fixed by `ruff check --fix`. Some errors will need to be fixed manually though.
- `black .` (notice the dot). This simply applies black formatting to the whole repo.
- `isort --profile black .` (notice the dot). This applies import sorting. You almost never need to run black/isort manually, pre-commit will take care of that.

#### MyPy

Python 3.5 introduced type annotations ([PEP 484](https://peps.python.org/pep-0484/)) and [MyPy](https://mypy.readthedocs.io/en/stable/) is the official type checker. You can add a VSCode/Cursor extension that will highlight most of the type issues. [Getting Started](https://mypy.readthedocs.io/en/stable/getting_started.html#) page in official docs is a good introduction to MyPy. Basically you specify the types for function arguments/outputs, which makes code easier to understand and helps catch issues with code early on.

### Commit Prefixes

- FEAT: - new features
- FIX: bug and other fixes on existing code
- DOCS: - any changes to readme, comments, docstrings etc
- STYLE: any style improvements, like refactoring, adding type hints, better namings etc
- DEV: intermediate commits, ideally these are on non-main branch. Basically every code after FEAT/FIX/DOCS/STYLE should be assumed to work correctly (to the best of your knowledge). DEV Is commits you do when working towards something that will become a FEAT/FIX on a main

## Acknowledgement

We made slight modifications to the project workflow, MkDocs setup, developer guides, and logging practices introduced by [Anton Morgunov](https://github.com/anmorgunov).
