"""
The MDForce Class
"""

import numpy as np

class MDForce:
    """ MDForce parent class with common methods and properties """
    @property
    def params(self):
        return self._params.copy()

    def __init__(self):
        self._params = dict()

    def set_params(self, **kwargs):
        self._params.update(kwargs)

    def get_params(self):
        return self._params.copy()

    def print_params(self):
        print(self._params)

    def compute_force(self, coords):
        raise NotImplementedError('compute_force should be implemented in child class')
