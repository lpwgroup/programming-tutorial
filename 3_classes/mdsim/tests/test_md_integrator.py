import numpy as np
from mdsim.integrator import MDIntegrator

def test_init():
    integrator = MDIntegrator()
    assert isinstance(integrator, MDIntegrator)