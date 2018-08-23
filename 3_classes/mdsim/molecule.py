"""
The Molecule Class
"""

import numpy as np

class Molecule:
    ELEM_MASS = {'H': 1.0, 'He': 2.0}

    @property
    def noa(self):
        return len(self.elems)

    @property
    def masses(self):
        return np.array([self.ELEM_MASS[e] for e in self.elems], dtype=float)

    def __init__(self):
        # elements list
        self.elems = []
        # coords
        self.coords = []
        # force
        self.force = []

    def set_elems(self, elems):
        assert all(e in self.ELEM_MASS for e in elems), 'Element not recognized'
        self.elems = elems

    def set_coords(self, coords):
        """ Set coordinates of molecule """
        assert coords.shape == (self.noa, 3)
        self.coords = coords

    def set_force(self, force):
        """ Set force of molecule """
        assert force.shape == (self.noa, 3)
        self.force = force

    def reset_force(self):
        self.force = np.zeros((self.noa, 3), dtype=float)

    def create_cube(self, n=3, element='He'):
        """ Create a molecule as a cube of atoms.

        Parameters
        ----------
        n: integer
            number of atoms in each dimension of the cube
        element: string
            The element of all atoms in this molecule
        """
        self.set_elems([element] * n**3)
        self.set_coords(np.array([[x,y,z] for x in range(n) for y in range(n) for z in range(n)], dtype=float))

    def copy(self):
        """ Return a new copy of self """
        new_molecule = Molecule()
        new_molecule.elems = self.elems.copy()
        new_molecule.coords = self.coords.copy()
        new_molecule.force = self.force.copy()
        return new_molecule