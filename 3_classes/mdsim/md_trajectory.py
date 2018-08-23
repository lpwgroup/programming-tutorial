"""
The MDTrajectory Class
"""

import numpy as np
from mdsim.molecule import Molecule

class MDTrajectory:
    """
    MDTrajectory is created for reading and writing trajectory files, and
    provides some analyze methods on the trajectory.
    """
    @property
    def xyz(self):
        return np.array(self._xyz, dtype=float)

    @property
    def noa(self):
        return len(self.elems)

    @property
    def length(self):
        return len(self._xyz)

    def __len__(self):
        return len(self._xyz)

    def __init__(self, elems):
        self._xyz = []
        self.elems = elems

    def add_frame(self, molecule):
        """ Add a frame in trajectory

        Parameters
        ----------
        molecule: Molecule object
            The molecule to be added as a frame
        """
        assert isinstance(molecule, Molecule)
        assert np.array_equal(self.elems, molecule.elems), f'{self.elems} != {molecule.elems}'
        self._xyz.append(molecule.coords.copy())

    def load_xyz(self, filename):
        """ Read the xyz file as a trajectory """
        # load all lines at once
        with open(filename) as f:
            lines = f.readlines()
        n_lines = len(lines)
        self._xyz = []
        # read the number of atoms from first line
        noa = int(lines[0])
        lines_per_frame = noa + 2
        # read the frames
        for i_frame in range(int(n_lines / lines_per_frame)):
            frame_lines = lines[i_frame*lines_per_frame:(i_frame+1)*lines_per_frame]
            frame_coords = self.parse_xyz_frame(frame_lines)
            self._xyz.append(frame_coords)

    def parse_xyz_frame(self, frame_lines):
        """ Parse a list of xyz lines into a numpy array of coordinates """
        coords = []
        for line in frame_lines[2:]:
            coords.append(np.array(line.split()[1:], dtype=float))
        return np.array(coords, dtype=float)

    def save_xyz(self, filename):
        """ Write the trajectory in xyz format """
        with open(filename, 'w') as outfile:
            for i, frame_xyz in enumerate(self._xyz):
                outfile.write(f'{self.noa}\nFrame{i:10d}\n')
                for e, c in zip(self.elems, frame_xyz):
                    outfile.write(f'{e} {c[0]:10.7f} {c[1]:10.7f} {c[2]:10.7f}\n')

    def find_break_frame(self, thresh=0.1):
        """ Find the frame in trajectory where the structure changes from initial
        Parameters
        ----------
        thresh: float
            The threshold for detecting the geometry change

        Returns
        -------
        break_frame: int
            The index of the first frame that the geometry changes
        """
        xyz = self.xyz
        # compute max abs displacements along the trajectory
        max_move = np.abs(xyz - xyz[0]).max(axis=(1,2))
        # find the first frame that displacement is larger than thresh
        break_frame = (max_move > thresh).argmax()
        return break_frame
