# Python MD Program written in functions

It is highly suggested to put codes into functions in any case when the structure is not "flat".

There are many benefits for doing this:

- Functions calling is a natural way to think about the work flow.

- Intuitive function names, input and output parameters are very helpful for understanding the code.

- Many Python features are enabled with functions. (import, testing, profiling, decorator etc.)

## Comparing to the flat script

- Improved LJ force calculation enabled with a new function.

- Analysis function imported from separate file.

- cProfile shows how long each function takes.

    `python -m cProfile -s time md_run.py`

- Argument parser implemented.

- Tests with `pytest -v`

## Files in this folder

1. `md_run.py`

    The main python code implemented the molecular dynamics (MD) program.

    Notes:

    - main() function

    - Function doc string in Numpy format

2. `md_analysis.py`

    The python code that loads a trajectory and run analysis.

    Notes:

    -  Functions in it can be imported to run in other codes

3. `test_md.py`

    The file contain tests of many functions implemented in other codes.
