import numpy as np

__all__ = ['Hernquist']

class Hernquist:

    def __init__(self, unit, Rs, rho0=None, sigma0=None):
        if unit == 'arcsec':
            if not ((rho0 == None) and (sigma0 != None)):
                raise ValueError('When unit is arcsec, input only Rs and sigma0 in angular units!')
            rho0 = self.sigma02rho(sigma0=sigma0, Rs=Rs)
        elif unit == 'pix':
            if not ((rho0 != 0) and (sigma0 == None)):
                raise ValueError('When unit is pix, input only Rs and rho0 in pixel unit!')
        else:
            raise ValueError('Available units: [arcsec, pix].')
        self.rho0 = rho0
        self.Rs = Rs

    @staticmethod
    def sigma02rho(sigma0, Rs):
        """convert sigma0 to 3d density normalisation

        Args:
            sigma0 (float): projected surface density [Sigma_crit]
            Rs (float): scale radius [arcsec]

        Returns:
            float: 3d density normalisation [Sigma_crit arcsec^-1]
        """
        return sigma0 / Rs

    def Density_in_project(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**3)

