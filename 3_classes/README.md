# Python MD Simulation Package

Objective-oriented programming is suitable for building a "package" or "library".

The development is around classes. The class instance is a wrapped "state" with methods.

This type fully adopts the modularization, give maximum extendability and greatly simplify the maintenance.

Useful links:

- SOLID principle: https://deviq.com/solid/

- Design patterns: https://www.toptal.com/python/python-design-patterns

- MolSSI best practices: https://molssi.org/education/best-practices/

## Files in this folder

1. `mdsim/`

    The package root folder.

    `mdsim/md_simulation.py`: The `MDSimulation` class

    `mdsim/force/`: The `force` package

    `mdsim/integrator/`: The `integrator` package

    `mdsim/md_trajectory.py`: The `MDTrajectory` class

    `mdsim/molecule.py`: The `Molecule` class

    `mdsim/__init__.py`: required for import

2. `examples/run_mdsim.py`

    The script running a simulation utilizing the `mdsim` package.

3. `setup.py`

    The setup script that enables installation by `python setup.py install`


## Requirements

- Python 3.6 or above
- Numpy
- pytest
- setuptools
