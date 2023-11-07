import numpy as np
import astropy.constants as const
import astropy.units as u

__all__ = ['SIS_truncated_angular']

rad2arcsec = 360 * 3600 / 2 / np.pi

class SIS_truncated_angular(object):

    def __init__(self, sigma_v, lens_cosmo, rc = None):
        """The "truncated" SIS profile class

        Args:
            sigma_v (_type_): velocity dispersion [km/s]
            rc (_type_): truncation ardius [arcsec]
            lens_cosmo (_type_): LensCosmo instance
        """
        if rc == None:
            rc = 1000 * lens_cosmo.sis_sigma_v2theta_E(sigma_v)

        self.sigma_v = sigma_v
        self.rc = rc
        self.lens_cosmo = lens_cosmo
    
    def Density_3d_spherical(self, r):
        """3d density

        Args:
            r (_type_): radius [arcsec]

        Returns:
            _type_: _description_
        """
        return (((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec**2 / (r**2 * (1 + r**2/self.rc**2))

    
    def Density_3d_triaxial(self, x, y, z, zeta, xi, get_effective_radius = False):
        """3d density

        Args:
            x (_type_): x coordinate, aligned with the major axis
            y (_type_): y coordinate, aligned with the intermediate axis
            z (_type_): z coordinate, aligned with the minor axis
            zeta (_type_): intrinsic axis ratio b/a
            xi (_type_): intrinsic axis ratio c/a
            get_effective_radius (bool, optional): whether to return effective radius. Defaults to False.

        Returns:
            _type_: _description_
        """
        av_sq = (zeta * xi)**(2/3) * (x**2 + y**2/zeta**2 + z**2/xi**2)
        density = (((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec**2 / (av_sq * (1 + av_sq/self.rc**2))

        if get_effective_radius:
            return density, np.sqrt(av_sq)
        else:
            return density


    def Project_integrand(self, z, as_sq):
        """_summary_

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            float: return intrinsic density as integrand for 1d/2d projections (projected density is in [M_sun/kpc^2]) 
        """
        return (((self.sigma_v * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * rad2arcsec * (self.lens_cosmo.Dd * 1e3) / ((z**2 + as_sq) * (1 + (z**2 + as_sq)/self.rc**2))
        