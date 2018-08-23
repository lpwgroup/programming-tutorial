import pytest
import numpy as np
from mdsim.molecule import Molecule
from mdsim.integrator import MDIntegrator

def test_init():
    integrator = MDIntegrator()
    assert isinstance(integrator, MDIntegrator)

def test_integrate():
    integrator = MDIntegrator()
    molecule = Molecule()
    # integrate() is not implemented in parent class
    with pytest.raises(NotImplementedError):
        integrator.integrate(molecule)