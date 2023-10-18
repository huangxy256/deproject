import numpy as np
import astropy.constants as const
import astropy.units as u

__all__ = ['SIS_truncated_angular']

rad2arcsec = 360 * 3600 / 2 / np.pi

class SIS_truncated_angular(object):

    def __init__(self, sigma, rc, lens_cosmo):

        self.sigma = sigma
        self.rc = rc
        self.lens_cosmo = lens_cosmo
    
    def Density_3d_spherical(self, r):
        """Compute the spherical density for truncated SIS profile; the density follows: \rho \propto 1/(r**2 * (1 + r**2/rc**2))

        Args:
            r (_type_): _description_
            r_unit (_type_): _description_
            lens_cosmo (_type_, optional): _description_. Defaults to None.

        Raises:
            ValueError: _description_
            ValueError: _description_

        Returns:
            _type_: _description_
        """
            
        return (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec**2 / (r**2 * (1 + r**2/self.rc**2))

    
    def Density_3d_triaxial(self, x, y, z, zeta, xi, get_effective_radius = False):

        av_sq = (zeta * xi)**(2/3) * (x**2 + y**2/zeta**2 + z**2/xi**2)
        density = (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec**2 / (av_sq * (1 + av_sq/self.rc**2))

        if get_effective_radius:
            return density, np.sqrt(av_sq)
        else:
            return density


    def Density_in_project(self, z, as_sq):
        """_summary_

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            float: return intrinsic density as integrand for 1d/2d projections (projected density is in [M_sun/kpc^2]) 
        """
        return (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec * (self.lens_cosmo.Dd * 1e3) / ((z**2 + as_sq) * (1 + (z**2 + as_sq)/self.rc**2))
        