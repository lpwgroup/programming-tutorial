import pytest
import numpy as np
from mdsim.force import MDForce, LJForce

def test_init():
    f = LJForce()
    assert isinstance(f, LJForce)
    # LJForce instance should also be an MDForce instance
    assert isinstance(f, MDForce)

def test_set_params():
    f = LJForce()
    f.set_params(sigma=0.9, epsilon=1.0)
    assert f.params['sigma'] == 0.9
    assert f.params['epsilon'] == 1.0
    # directly overwrite f.params should not be allowed
    with pytest.raises(AttributeError):
        f.params = {'a': 1}

def test_get_params():
    f = LJForce()
    f.set_params(sigma=0.9, epsilon=1.0)
    assert f.get_params() == f.params

def test_compute_force_ref():
    f = LJForce()
    f.set_params(sigma=0.9, epsilon=1.0)
    coords = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    # check if compute_force_ref gives correct result
    ref_force = f.compute_force_ref(coords)
    f0 = np.array([
        [ 0.80203375,  0.80203375,  0.80203375],
        [-0.69125042, -0.69125042,  0.5804671 ],
        [-0.69125042,  0.5804671 , -0.69125042],
        [ 0.5804671 , -0.69125042, -0.69125042]
    ])
    assert np.allclose(ref_force, f0)

def test_compute_force():
    """ test LJForce compute_force against compute_force_ref """
    f = LJForce()
    f.set_params(sigma=0.9, epsilon=1.0)
    coords = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    for _ in range(10):
        rand_coords = coords + 0.5 * np.random.random(coords.shape)
        force = f.compute_force(rand_coords)
        ref_force = f.compute_force_ref(rand_coords)
        assert np.allclose(force, ref_force)
