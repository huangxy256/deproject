import numpy as np

from scipy.integrate import quad

__all__ = ['Hernquist']

class Hernquist:

    def __init__(self, rho0, Rs):
        self.rho0 = rho0
        self.Rs = Rs

    def Density_in_project(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**3)



