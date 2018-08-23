import numpy as np
from mdsim.force.md_force import MDForce

class LJForce(MDForce):
    """
    LJForce class implements the Lennard-Jones force

    Reference
    ---------
    The LJ force takes the formula:
        Fij = (-12 x sigma^12 / rij^13 + 6 x sigma^6 / rij^7) * 4 * epsilon * [rij]/rij
    """
    def __init__(self):
        super().__init__()
        # set default sigma and epsilon parameter
        self._params = {'sigma': 1.0, 'epsilon': 1.0}

    def compute_force(self, coords):
        """
        Compute the Lennard-Jones force vector on current geometry

        Parameters
        ----------
        coords: numpy.ndarray of shape (Natoms, 3)
            Numpy array of atomic coordinates

        Returns
        -------
        force: numpy.ndarray of shape (Natoms, 3)
            Numpy array of gradients on each atom
        """
        sigma = self._params['sigma']
        epsilon = self._params['epsilon']
        # compute the distance between each atom pairs
        c_diff = coords[:,np.newaxis,:] - coords[np.newaxis,:,:]
        r2_mat = np.sum(np.square(c_diff), axis=-1)
        # prepare values for the LJ force formula
        s6 = sigma**6
        r2_mat2 = np.square(r2_mat)
        np.fill_diagonal(r2_mat2, 1.0) # prevent 1/0 error in equation
        r2_matn4 = 1.0 / np.square(r2_mat2)
        r2_matn7 = np.square(r2_matn4) * r2_mat
        # compute the magnitude of the gradients
        f_lj_mat = (-12.0*r2_matn7 * s6 + 6.0*r2_matn4) * 4 * epsilon * s6
        # contract with coordinates dR to get gradient vectors
        return np.einsum('ijk,ij->ik', c_diff, f_lj_mat)

    def compute_force_ref(self, coords):
        """
        Compute the Lennard-Jones force vector on current geometry

        Parameters
        ----------
        coords: numpy.ndarray of shape (Natoms, 3)
            Numpy array of atomic coordinates

        Returns
        -------
        force: numpy.ndarray of shape (Natoms, 3)
            Numpy array of gradients on each atom
        """

        noa = len(coords)
        s6 = self._params['sigma']**6
        epsilon = self._params['epsilon']
        forces = np.zeros((noa,3), dtype=float)
        for i in range(noa):
            for j in range(i+1,noa):
                dc = coords[i] - coords[j]
                r2 = dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]
                f = (-12 / r2**7 * s6 + 6 / r2**4) * 4 * epsilon * s6
                forces[i] += f * dc
                forces[j] -= f * dc
        return forces