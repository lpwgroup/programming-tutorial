import numpy as np
from mdsim.molecule import Molecule
from mdsim.md_trajectory import MDTrajectory

def test_init():
    elems = ['He', 'He']
    traj = MDTrajectory(elems)
    assert isinstance(traj, MDTrajectory)
    assert traj.elems == elems
    assert traj.noa == 2
    assert len(traj) == traj.length == 0

def test_add_frame():
    elems = ['He', 'He']
    coords = np.array([[0,0,0], [0,0,1]])
    molecule = Molecule()
    molecule.set_elems(elems)
    molecule.set_coords(coords)
    traj = MDTrajectory(elems)
    traj.add_frame(molecule)
    assert len(traj) == 1

def test_fine_break_frame():
    elems = ['He', 'He']
    coords = np.array([[0,0,0], [0,0,1]])
    molecule = Molecule()
    molecule.set_elems(elems)
    molecule.set_coords(coords)
    traj = MDTrajectory(elems)
    # add 5 of the same frames
    for _ in range(5):
        traj.add_frame(molecule)
    # add a different 6th frame
    molecule1 = molecule.copy()
    molecule1.set_coords(coords + 1)
    traj.add_frame(molecule1)
    break_frame = traj.find_break_frame(thresh=0.1)
    assert break_frame == 5