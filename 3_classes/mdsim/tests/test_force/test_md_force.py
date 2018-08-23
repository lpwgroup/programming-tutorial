import pytest
import numpy as np
from mdsim.force import MDForce

def test_init():
    f = MDForce()
    assert isinstance(f, MDForce)

def test_set_params():
    f = MDForce()
    f.set_params(a=1, b=2)
    assert f.params['a'] == 1
    assert f.params['b'] == 2
    # directly overwrite f.params should not be allowed
    with pytest.raises(AttributeError):
        f.params = {'a': 1}

def test_get_params():
    f = MDForce()
    f.set_params(a=1, b=2)
    assert f.get_params() == f.params

def test_compute_force():
    f = MDForce()
    f.set_params(sigma=0.9, epsilon=1.0)
    coords = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    # compute force is not implemented in parent class
    with pytest.raises(NotImplementedError):
        f.compute_force(coords)
