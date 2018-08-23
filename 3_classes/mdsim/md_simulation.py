"""
The MDSimulation Class
"""

import numpy as np
from mdsim.molecule import Molecule
from mdsim.force import MDForce
from mdsim.integrator import MDIntegrator
from mdsim.md_trajectory import MDTrajectory

class MDSimulation:
    """ The MDSimulation class wrapps other modules and carry out simulations """
    @property
    def trajectory(self):
        return self._traj

    def __init__(self, molecule, ff, integrator, interval=100, verbose=False):
        assert isinstance(molecule, Molecule)
        assert isinstance(ff, MDForce)
        assert isinstance(integrator, MDIntegrator)
        self.molecule = molecule
        self.ff = ff
        self.integrator = integrator
        self.interval = interval
        self.verbose = verbose
        self._traj = MDTrajectory(molecule.elems)
        self._step = 0

    def step(self, n_steps):
        """ Run simulation for n_steps """
        for _ in range(n_steps):
            force = self.ff.compute_force(self.molecule.coords)
            self.molecule.set_force(force)
            self.molecule = self.integrator.integrate(self.molecule)
            if self._step % self.interval == 0:
                self._traj.add_frame(self.molecule)
                # print verbose information
                if self.verbose:
                    print(f"step {self._step:15d}")
            # increment step count
            self._step += 1
