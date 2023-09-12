import numpy as np

__all__ = ['NFW']

class NFW(object):

    def __init__(self, Rs, unit, rho0=None, alpha_Rs=None):
        if unit == 'arcsec':
            if not ((rho0 == None) and (alpha_Rs != None)):
                raise ValueError('When unit is arcsec, input only Rs and alpha_Rs in arcsec!')
            rho0 = self.alpha2rho0(alpha_Rs=alpha_Rs, Rs=Rs)
        elif unit == 'pix':
            if not ((rho0 != 0) and (alpha_Rs == None)):
                raise ValueError('When unit is pix, input only Rs and rho0 in pixel unit!')
        else:
            raise ValueError('Available units: [arcsec, pix].')
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

    def Density_in_project(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**2)