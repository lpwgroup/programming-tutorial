#!/usr/bin/env python

#####################
# Simple MD program #
#####################

import time
import numpy as np

def create_molecule(n=3, element='He'):
    """
    Create a molecule as atoms in a cube.

    Parameters
    ----------
    n: integer
        number of atoms in each dimension of the cube
    element: string
        The element of all atoms in this molecule

    Returns
    -------
    coords: numpy.ndarray of shape (n**3, 3)
        Numpy array of atomic coordinates
    elems: list of strings
        List of elements for all atoms
    """
    coords = np.array([[x,y,z] for x in range(n) for y in range(n) for z in range(n)], dtype=float)
    elems = [element] * len(coords)
    return coords, elems

def ref_LJ_force(coords, epsilon=1.0, sigma=1.0):
    """ Compute the LJ force for the molecule in current geometry.

    Parameters
    ----------
    coords: Numpy.ndarray of shape (Natoms, 3)
        Numpy array of atomic coordinates
    epsilon: float
        epsilon parameter in LJ force formula
    sigma: float
        sigma parameter in LJ force formula

    Returns
    -------
    force: Numpy.ndarray of shape (Natoms, 3)
        Numpy array of gradients on each atom

    Reference
    ---------
    The LJ force takes the formula:
        Fij = (-12 x sigma^12 / rij^14 + 6 x sigma^6 / rij^8) * 4 * epsilon
    """
    noa = len(coords)
    s6 = sigma**6
    forces = np.zeros((noa,3))
    for i in range(noa):
        for j in range(i+1,noa):
            dc = coords[i] - coords[j]
            r2 = dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]
            f = (-12 / r2**7 * s6 + 6 / r2**4) * 4 * epsilon * s6
            forces[i] += f * dc
            forces[j] -= f * dc
    return forces

def numpy_LJ_force(coords, epsilon=1.0, sigma=1.0):
    """ Compute the LJ force for the molecule in current geometry.

    Parameters
    ----------
    coords: Numpy.ndarray of shape (Natoms, 3)
        Numpy array of atomic coordinates
    epsilon: float
        epsilon parameter in LJ force formula
    sigma: float
        sigma parameter in LJ force formula

    Returns
    -------
    force: Numpy.ndarray of shape (Natoms, 3)
        Numpy array of gradients on each atom

    Reference
    ---------
    The LJ force takes the formula:
        Fij = (-12 x sigma^12 / rij^14 + 6 x sigma^6 / rij^8) * 4 * epsilon
    """
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

ELEM_MASS = {'He': 2}

def verlet_intergrate(forcefunc, coords, elems, nsteps=10000, dt=0.001, outfile='traj.xyz', verbose=True):
    """ Run MD simulation with Verlet intergration.

    Parameters
    ----------
    forcefunc: function for computing force
        This function should take coords, sigma and epsilon as input, then return a numpy array in same shape as coords
    coords: numpy.ndarray of shape (Natoms, 3)
        Initial coordinates of the molecule
    elems: list of strings
        List of elements for all atoms
    nsteps: int
        Number to total steps in MD simulation
    step_length: float
        Time length of each MD step
    outfile: str
        The name of the output file
    verbose: bool
        Flag for printing verbose information

    Returns
    -------
    traj: numpy.ndarray of shape (Nframes, Natoms, 3)
        The trajectory of the simulation
    """
    # traj contains the list of frames
    traj = []
    prev_coords = coords.copy()
    masses = np.array([ELEM_MASS[e] for e in elems])
    with open(outfile,'w') as outputfile:
        for i_frame in range(nsteps):
            # compute force
            force = forcefunc(coords, sigma=0.9, epsilon=20.0)
            dr = -force / masses[:, np.newaxis] * dt**2
            #update coords
            tmp = coords.copy()
            coords = 2 * coords + dr - prev_coords
            prev_coords = tmp
            if i_frame % 100 == 0:
                # append the coords to traj list
                traj.append(coords.copy())
                # write the coordinate to file in xyz format
                write_xyz_frame(coords, elems, i_frame, outputfile)
                if verbose:
                    print(f'{i_frame:6d} steps finished', end='\r', flush=True)
    return np.array(traj, dtype=float)

def write_xyz_frame(coords, elems, frame_number, outputfile):
    """ Write the current frame as xyz format into an opened file handle """
    outputfile.write(f'{len(coords)}\nframe {frame_number}\n')
    for c,e in zip(coords, elems):
        outputfile.write(f'{e} {c[0]:10.7f} {c[1]:10.7f} {c[2]:10.7f}\n')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Simple MD program', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--atoms_each_edge', type=int, default=3, help="Number of atoms in each edge of cube.")
    parser.add_argument('-e', '--element', type=str, default='He', help="Element of all atoms.")
    parser.add_argument('-n', '--num_steps', type=int, default=10000, help="Number of intergration steps.")
    parser.add_argument('-f', '--force_function', type=str, choices=['ref', 'numpy'], default='ref', help="Function to compute force")
    parser.add_argument('-o', '--outfile', type=str, default='traj.xyz', help='Name of output trajectory file')
    parser.add_argument('--analysis', action='store_true', help='Flag to enable analysis following simulation')
    args = parser.parse_args()

    # step 1: create molecule
    coords, elems = create_molecule(args.atoms_each_edge)

    # step 2: select force function
    ffdict = {
        'ref': ref_LJ_force,
        'numpy': numpy_LJ_force,
    }
    ffunc = ffdict[args.force_function]

    # step 3: run intergrator
    traj = verlet_intergrate(ffunc, coords, elems, nsteps=args.num_steps, outfile=args.outfile, verbose=True)

    # step 4: analyze trajectory
    if args.analysis:
        from md_analysis import find_break_frame
        break_frame = find_break_frame(traj) * 100
        print(f"Found cube destructs at frame ~ {break_frame}")

if __name__ == '__main__':
    main()