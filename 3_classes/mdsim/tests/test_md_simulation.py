from mdsim.molecule import Molecule
from mdsim.force import LJForce
from mdsim.integrator import VerletIntegrator
from mdsim.md_trajectory import MDTrajectory
from mdsim.md_simulation import MDSimulation

def test_init():
    molecule = Molecule()
    molecule.create_cube(n=3)
    ff = LJForce()
    ff.set_params(sigma=0.9, epsilon=20.0)
    integrator = VerletIntegrator(t_step=0.001)
    simulation = MDSimulation(molecule, ff, integrator, interval=100)
    assert isinstance(simulation, MDSimulation)
    assert isinstance(simulation.trajectory, MDTrajectory)

def test_step():
    # create simulation
    molecule = Molecule()
    molecule.create_cube(n=3)
    ff = LJForce()
    ff.set_params(sigma=0.9, epsilon=20.0)
    integrator = VerletIntegrator(t_step=0.001)
    simulation = MDSimulation(molecule, ff, integrator, interval=100)
    # test running simulation for 100 steps
    simulation.step(100)
    assert simulation.current_step == 100