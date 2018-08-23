"""
The MDIntegrator Class
"""

class MDIntegrator:
    """ The MDIntegrator class with common methods and properties """
    def __init__(self, t_step=0.001):
        self.t_step = t_step

    def integrate(self, molecule):
        raise NotImplementedError("integrate() method should be implemented in child class")
