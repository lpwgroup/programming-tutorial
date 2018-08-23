#!/usr/bin/env python

import numpy as np

def load_xyz_traj(filename):
    """ Read the xyz file as a trajectory """
    # load all lines at once
    with open(filename) as f:
        lines = f.readlines()
    n_lines = len(lines)
    traj = []
    # read the number of atoms from first line
    noa = int(lines[0])
    lines_per_frame = noa + 2
    # read the frames
    for i_frame in range(int(n_lines / lines_per_frame)):
        frame_lines = lines[i_frame*lines_per_frame:(i_frame+1)*lines_per_frame]
        frame_coords = parse_xyz_frame(frame_lines)
        traj.append(frame_coords)
    return np.array(traj, dtype=float)

def parse_xyz_frame(lines):
    """ read the traj coordinates from lines for one frame """
    coords = []
    for line in lines[2:]:
        coords.append(np.array(line.split()[1:], dtype=float))
    return np.array(coords)

def find_break_frame(traj, thresh=0.1):
    """ Find the frame in trajectory where the structure changes from initial
    Parameters
    ----------
    traj: numpy.ndarray of shape (Nframes x Natoms x 3)
        Coordinates of atoms in a md trajectory
    thresh: float
        The threshold for detecting the geometry change

    Returns
    -------
    break_frame: int
        The index of the first frame that the geometry changes
    """
    # compute max abs displacements along the trajectory
    max_move = np.abs(traj - traj[0]).max(axis=(1,2))
    # find the first frame that displacement is larger than thresh
    break_frame = (max_move > thresh).argmax()
    return break_frame

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Analyze MD trajectory', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('trajfile', help="Input trajectory file.")
    parser.add_argument('-t', '--thresh', type=float, default=0.1, help="Threshold for detecting the geometry change")
    args = parser.parse_args()

    traj = load_xyz_traj(args.trajfile)
    break_frame = find_break_frame(traj)

    print(f"Found cube destructs at frame ~ {break_frame}")

if __name__ == '__main__':
    main()