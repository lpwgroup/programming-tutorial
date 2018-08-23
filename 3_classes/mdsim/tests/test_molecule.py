import numpy as np
from mdsim.molecule import Molecule

def test_init():
    m = Molecule()
    assert m.noa == 0

def test_set_elems():
    m = Molecule()
    m.set_elems(['He', 'He'])
    assert m.noa == 2
    assert np.array_equal(m.masses, [2.0,2.0])

def test_set_coords():
    m = Molecule()
    m.set_elems(['He', 'He'])
    rand_coords = np.random.random((2,3))
    m.set_coords(rand_coords)
    assert np.array_equal(m.coords, rand_coords)

def test_set_force():
    m = Molecule()
    m.set_elems(['He', 'He'])
    force = np.random.random((2,3))
    m.set_force(force)
    assert np.array_equal(m.force, force)

def test_reset_force():
    m = Molecule()
    m.set_elems(['He', 'He'])
    force = np.random.random((2,3))
    m.set_force(force)
    assert np.array_equal(m.force, force)
    m.reset_force()
    assert np.array_equal(m.force, np.zeros((2,3)))

def test_create_cube():
    for na in range(1,5):
        m = Molecule()
        m.create_cube(n=na, element='He')
        assert m.coords.shape == (na**3, 3)
        assert m.elems == ['He'] * na**3

def test_copy():
    m = Molecule()
    m.create_cube(n=3, element='He')
    force = np.random.random((3**3,3))
    m.set_force(force)
    m1 = m.copy()
    assert np.array_equal(m.elems, m1.elems)
    assert np.array_equal(m.masses, m1.masses)
    assert np.array_equal(m.coords, m1.coords)
    assert np.array_equal(m.force, m1.force)
    # check if modify one changes the other
    m1.elems[0] = 'H'
    m1.coords += 1
    m1.force += 1
    assert not np.array_equal(m.elems, m1.elems)
    assert not np.array_equal(m.masses, m1.masses)
    assert not np.array_equal(m.coords, m1.coords)
    assert not np.array_equal(m.force, m1.force)
