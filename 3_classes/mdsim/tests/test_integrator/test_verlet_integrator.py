import pytest
import numpy as np
from mdsim.molecule import Molecule
from mdsim.integrator import MDIntegrator, VerletIntegrator

def test_init():
    integrator = VerletIntegrator()
    assert isinstance(integrator, VerletIntegrator)
    assert isinstance(integrator, MDIntegrator)

def test_integrate():
    integrator = VerletIntegrator(t_step=0.001)
    molecule = Molecule()
    molecule.set_elems(['He', 'He', 'He', 'He'])
    coords = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    molecule.set_coords(coords)
    force = np.array([
        [ 0.80203375,  0.80203375,  0.80203375],
        [-0.69125042, -0.69125042,  0.5804671 ],
        [-0.69125042,  0.5804671 , -0.69125042],
        [ 0.5804671 , -0.69125042, -0.69125042]
    ])
    molecule.set_force(force)
    # test integrate for the first time without prev_coords
    molecule1 = integrator.integrate(molecule)
    ref_disp = -1 * force / molecule.masses[:, np.newaxis] * 0.001**2
    assert np.allclose(molecule1.coords, coords + ref_disp)
    assert np.array_equal(integrator.prev_coords, coords)
    # test integrate for the second time with prev_coords
    molecule1.set_force(force)
    molecule2 = integrator.integrate(molecule1)
    ref_disp = -1 * force / molecule.masses[:, np.newaxis] * 0.001**2
    ref_disp += molecule1.coords - molecule.coords
    assert np.allclose(molecule2.coords, molecule1.coords + ref_disp)