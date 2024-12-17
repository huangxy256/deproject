import numpy as np

__all__ = ['Jaffe']

class Jaffe(object):

    def __init__(self, rs, rho0):
        """class to manage the Jaffe profile, which is rho(r) = rho0 * (r/rs)^(-2) * (1 + r/rs)^(-2)

        Args:
            rs (_type_): scale radius 
            rho0 (_type_): mass normalization parameter
        """
        self.rs = rs
        self.rho0 = rho0

    def Density_3d_spherical(self, r):
        """3d density 

        Args:
            r (_type_): radius [arcsec]

        Returns:
            _type_: _description_
        """
        return self.rho0 * (r / self.rs)**(-2) * (1 + (r / self.rs))**(-2)
    
    def Density_3d_triaxial(self, x, y, z, zeta, xi, get_effective_radius = False):
        """3d density

        Args:
            x (_type_): x coordinate, aligned with intrinsic major axis
            y (_type_): y coordinate, aligned with intrinsic intermediate axis
            z (_type_): z xoordinate, aligned with intrinsic minor axis
            zeta (_type_): intrinsic axis ratio b/a
            xi (_type_): intrinsic axis ratio c/a
            get_effective_radius (bool, optional): whether to return effective radius. Defaults to False.

        Returns:
            _type_: _description_
        """
        av = (zeta * xi)**(1/3) * np.sqrt(x**2 + y**2/zeta**2 + z**2/xi**2)
        density = self.rho0 * (av / self.rs)**(-2) * (1 + (av / self.rs))**(-2)

        if get_effective_radius:
            return density, av
        else:
            return density
        
    def Project_integrand(self, z, as_sq):
        """Integrand in projection integral, z to be integrated 

        Args:
            z (float): line of sight variable
            as_sq (float): square of radial variable in projection

        Returns:
            float: density in projection
        """
        return self.rho0 / ((np.sqrt(z**2 + as_sq)/ self.rs)**2 * (1 + (np.sqrt(z**2 + as_sq)/ self.rs))**2)