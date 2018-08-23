#!/usr/bin/env python

import md_run
import md_analysis
import numpy as np

def test_create_molecule():
    for na in range(1,5):
        coords, elems = md_run.create_molecule(n=na, element='Te')
        assert coords.shape == (na**3, 3)
        assert elems == ['Te'] * na**3

def test_ref_LJ_force():
    # test force of one atom to be zero
    coords, elems = md_run.create_molecule(n=1)
    force = md_run.ref_LJ_force(coords)
    assert np.array_equal(force, np.zeros((1,3)))
    # test force of 2x2x2 atom cube,
    coords, elems = md_run.create_molecule(n=2)
    force = md_run.ref_LJ_force(coords, epsilon=1.0, sigma=0.9)
    # test total force to be 0
    assert np.allclose(np.sum(force), 0)
    # test one force element
    assert np.allclose(force[0], [-0.73173237, -0.73173237, -0.73173237])

def test_numpy_LJ_force():
    for na in range(1,5):
        coords, elems = md_run.create_molecule(n=na)
        numpy_force = md_run.numpy_LJ_force(coords, epsilon=1.0, sigma=0.9)
        # test total force to be 0
        assert np.allclose(np.sum(numpy_force), 0)
        # test force against ref
        ref_force = md_run.ref_LJ_force(coords, epsilon=1.0, sigma=0.9)
        assert np.allclose(numpy_force, ref_force)
        # test additional random geo
        for _ in range(10):
            disp = np.random.random(coords.shape)
            rand_coords = coords + disp
            numpy_force = md_run.numpy_LJ_force(rand_coords, epsilon=1.0, sigma=0.9)
            ref_force = md_run.ref_LJ_force(rand_coords, epsilon=1.0, sigma=0.9)
            assert np.allclose(numpy_force, ref_force)

def test_parse_xyz_frame():
    lines = ['3', '', 'H 0 0 0', 'H 0 0 1', 'H 1 2 3']
    ref_coords = np.array([[0,0,0], [0,0,1], [1,2,3]])
    coords = md_analysis.parse_xyz_frame(lines)
    assert np.allclose(coords, ref_coords)

def test_find_break_frame():
    traj = np.zeros((100, 10, 3))
    traj[25:] = 1
    break_frame = md_analysis.find_break_frame(traj, thresh=0.1)
    assert break_frame == 25