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

    def Density_spherical(self, rad):
        """3d density

        Args:
            rad (_type_): radius [arcsec]

        Returns:
            _type_: density as a function of radius [Sigma_crit]
        """
        return self.rho0 / ((rad / self.Rs) * (1 + (rad / self.Rs)**2))

    def Density_triaxial(self, x, y, z, zeta, xi, get_effective_radius=False):
        """3d density

        Args:
            x (_type_): x coordinate, aligned with intrinsic major axis
            y (_type_): y coordinate, aligned with intrinsic intermediate axis
            z (_type_): z xoordinate, aligned with intrinsic minor axis
            zeta (_type_): intrinsic axis ratio b/a
            xi (_type_): intrinsic axis ratio c/a
            get_effective_radius (bool, optional): whether to return effective radius. Defaults to False.

        Returns:
            _type_: density as a function of effective radius [Sigma_crit]
        """
        av = (zeta * xi)**(1/3) * np.sqrt(x**2 + y**2/zeta**2 + z**2/xi**2)
        density = self.rho0 / ((av / self.Rs) * (1 + (av / self.Rs)**2))
        
        if get_effective_radius:
            return density, av
        else:
            return av

    def Project_integrand(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.Rs) * (1 + (np.sqrt(z**2 + as_sq)/ self.Rs))**2)