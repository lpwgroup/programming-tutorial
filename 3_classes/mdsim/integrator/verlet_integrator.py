"""
The VerletIntegrator Class
"""

import numpy as np
from mdsim.molecule import Molecule
from mdsim.integrator.md_integrator import MDIntegrator

class VerletIntegrator(MDIntegrator):
    """ VerletIntegrator that implements the verlet algorithm

    References
    ----------
    Verlet integration formula:
        X[i+1] = X[i] x 2 + dX - X[i-1]
    """
    def integrate(self, molecule):
        """
        integrate current coords and force to get new coordinates

        Parameters
        ----------
        molecule: mdsim.molecule.Molecule object
            molecule.coords, molecule.force, molecule.masses will be used.

        Returns
        -------
        new_molecule: mdsim.molecule.Molecule object
            The new molecule containing the new set of coordinates
        """
        assert isinstance(molecule, Molecule)
        # initialize prev_coords
        if not hasattr(self, 'prev_coords'):
            self.prev_coords = molecule.coords.copy()
        # compute step dX
        dx = -molecule.force / molecule.masses[:, np.newaxis] * self.t_step**2
        # create new molecule
        new_molecule = molecule.copy()
        new_molecule.reset_force()
        new_molecule.set_coords(molecule.coords * 2 + dx - self.prev_coords)
        # update self.prev_coords
        self.prev_coords = molecule.coords.copy()
        return new_molecule