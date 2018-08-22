#!/usr/bin/env python

#####################
# Simple MD program #
#####################

import time
import numpy as np

#------------------------------
#   I. Create molecular system
#------------------------------

# number of atoms in each dimension forming a cube
n = 3

# Numpy array of shape (n_atoms x 3) for the coordinates of each atom
coords = np.array([[x,y,z] for x in range(n) for y in range(n) for z in range(n)], dtype=float)

# List of elements for each atom
elems = ['He'] * len(coords)


#------------------------------
#  II. MD Simulation
#------------------------------

# number of md steps
n_steps = 10000
# number of atoms
noa = len(coords)
# LJ Force parameters
sigma = 0.9
epsilon = 1.0
# traj as a list of frames
traj = []
# open an xyz file for writing the coordinates
with open('traj.xyz', 'w') as outfile:
    # verlet algorithm uses two coords, save one here
    prev_coord = coords.copy()
    # loop over md steps
    for i_frame in range(n_steps):
        # calculate the LJ force
        # the LJ force take the formula:
        # Fij = (-12 x sigma^12 / rij^13 + 6 x sigma^6 / rij^7) * 4 * epsilon * [rij]/rij
        s6 = sigma**6
        forces = np.zeros((noa,3), dtype=float)
        for i in range(noa):
            for j in range(i+1,noa):
                dc = coords[i] - coords[j]
                r2 = dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]
                f = (-12 / r2**7 * s6 + 6 / r2**4) * 4 * epsilon * s6
                forces[i] += f * dc
                forces[j] -= f * dc
        # use the force to update the coordinates
        mass = 0.1 # mass of each atom
        dt = 0.001 # time step
        dr = -forces / mass * dt**2
        # save a copy of current coords
        tmp = coords.copy()
        # verlet update step
        coords = coords * 2 + dr - prev_coord
        # update prev_coord for next step
        prev_coord = tmp
        # save every 100 frames
        if i_frame % 100 == 0:
            # append the coords to traj list
            traj.append(coords.copy())
            # write the coordinate to file in xyz format
            outfile.write(f'{noa}\nframe {i_frame}\n')
            for e, c in zip(elems, coords):
                outfile.write(f'{e:7s} {c[0]:10.7f} {c[1]:10.7f} {c[2]:10.7f}\n')
            # print the running status
            print(f'{i_frame:6d} steps finished', end='\r', flush=True)
print(f"MD Simulation of {n_steps} steps finished.")

#------------------------------
# III. Analyze the trajectory
#------------------------------

# convert the traj list to a numpy array
traj = np.array(traj)
# find the first frame that the cube destructs
max_move = np.abs(traj - traj[0]).max(axis=(1,2))
first_frame = (max_move > 0.1).argmax() * 100
print(f"Found cube destructs at frame ~ {first_frame}")