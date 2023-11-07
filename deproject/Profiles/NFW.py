import numpy as np

__all__ = ['NFW']

class NFW(object):

    def __init__(self, alpha_Rs, Rs):
        """the NFW profile class

        Args:
            alpha_Rs (_type_): _description_
            Rs (_type_): _description_
        """
        rho0 = self.alpha2rho0(alpha_Rs=alpha_Rs, Rs=Rs)
        self.rho0 = rho0
        self.Rs = Rs

    @staticmethod
    def alpha2rho0(alpha_Rs, Rs):
        """convert alpha_Rs to rho0

        Args:
            alpha_Rs (float): deflection angle at Rs [arcsec]
            Rs (float): scale radius [arcsec]

        Returns:
            float: density normalisation [Sigma_crit/arcsec]
        """
        rho0 = alpha_Rs / (4. * Rs**2 * (1. + np.log(0.5)))
        return rho0

    def Density_spherical(self, R):
        return 0

    def Density_triaxial(self):
        return 0

    def Project_integrand(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**2)