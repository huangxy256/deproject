import numpy as np
import astropy.constants as const
import astropy.units as u

__all__ = ['SIS']

class SIS(object):

    def __init__(self, sigma, lens_cosmo = None):

        self.sigma = sigma
        if lens_cosmo != None:
            self.lens_cosmo = lens_cosmo
        else:
            pass
    
    def Density_3d(self, r, r_unit, lens_cosmo = None):
        """_summary_

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
        if r_unit == 'kpc':
            return ((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (r * u.kpc)**2).to(u.M_sun / u.kpc**3).value
        elif r_unit == 'arcsec':

            if lens_cosmo == None:
                raise ValueError('No lens_cosmo instance provided. Please pass LensCosmo instance as a variable to this function.')
            else:
                pass
            
            return (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * (360 * 3600 / 2 / np.pi)**2 / r**2
        else:
            raise ValueError("Supported r_unit: ['kpc', 'arcsec'].")

    def _Density_in_project(self, z, as_sq):
        """Integrand of density in LoS integration [M_sun/kpc^3]

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            _type_: _description_
        """
        return ((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G * (z**2 + as_sq) * u.kpc**2)).to(u.M_sun / u.kpc**3).value

    def Density_in_project(self, z, as_sq):
        """_summary_

        Args:
            z (_type_): _description_
            as_sq (_type_): _description_

        Returns:
            float: return intrinsic density as integrand for 1d/2d projections (projected density is in [M_sun/kpc^2]) 
        """
        return (((self.sigma * u.km / u.s)**2 / (2 * np.pi * const.G) / (self.lens_cosmo.Dd * u.Mpc)**2).to(u.M_sun / u.kpc**3)).value * (360 * 3600 / 2 / np.pi) * (self.lens_cosmo.Dd * 1e3) / (z**2 + as_sq)
        