import numpy as np
import astropy.constants as const
import astropy.units as u

__all__ = ['SIS_truncated_physical']

class SIS_truncated_physical(object):

    def __init__(self, sigma_v, rc):
        """the "truncated" SIS profile class

        Args:
            sigma_v (_type_): _description_
            rc (_type_): _description_
        """
        self.sigma_v = sigma_v
        self.rc = rc
    
    def Density_3d_spherical(self, r):
        """Compute the spherical density for truncated SIS profile; the density follows: \rho \propto 1/(r**2 * (1 + r**2/rc**2))

        Args:
            r (_type_): _description_

        Raises:
            ValueError: _description_
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        return ((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G) / ((r * u.kpc)**2 * (1 + r**2/self.rc**2))).to(u.M_sun / u.kpc**3).value


    def Density_3d_triaxial(self, x, y, z, zeta, xi, get_effective_radius = False):
        """Compute the spherical density for truncated SIS profile; the density follows: \rho \propto 1/(r**2 * (1 + r**2/rc**2))

        Args:
            r (_type_): _description_

        Raises:
            ValueError: _description_
            ValueError: _description_

        Returns:
            _type_: _description_
        """

        av_sq = (zeta * xi)**(2/3) * (x**2 + y**2/zeta**2 + z**2/xi**2)
        density = ((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G) / (av_sq * u.kpc**2 * (1 + av_sq/self.rc**2))).to(u.M_sun / u.kpc**3).value
        if get_effective_radius:
            return density, np.sqrt(av_sq)
        else:
            return density
        

    def Project_integrand(self, z, as_sq):
        """Integrand of density in LoS integration [M_sun/kpc^3]

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            _type_: _description_
        """
        return ((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G * (z**2 + as_sq) * u.kpc**2 * (1 + (z**2 + as_sq)/self.rc**2))).to(u.M_sun / u.kpc**3).value


        