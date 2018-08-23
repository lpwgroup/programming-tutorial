#!/usr/bin/env python

from mdsim.molecule import Molecule
from mdsim.force import LJForce
from mdsim.integrator import VerletIntegrator
from mdsim.md_simulation import MDSimulation

# create initial molecule as a 3x3x3 cube
molecule = Molecule()
molecule.create_cube(n=3)

# create force field with LJ parameters
ff = LJForce()
ff.set_params(sigma=0.9, epsilon=20.0)

# create integrator
integrator = VerletIntegrator(t_step=0.001)

# create simulation
simulation = MDSimulation(molecule, ff, integrator, interval=100, verbose=True)

# run simulation
simulation.step(10000)

# save simulation trajectory
simulation.trajectory.save_xyz('traj.xyz')

# analyze simulation trajectory
break_frame = simulation.trajectory.find_break_frame(thresh=0.1)
print(f"Found cube destructs at frame ~ {break_frame}")