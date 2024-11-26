from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np

__all__ = ['Projection']

class Projection(object):
    def __init__(self, zeta, xi, theta, phi):
        """class to manage the projection of a triaxial halo

        Args:
            zeta (_type_): intrinsic axis ratio between intermediate and major axis, zeta < 1
            xi (_type_): intrinsic axis ratio between minor and major axis, xi < 1
            theta (_type_): polar angle in a spherical coordinate whose origin is at the halo's center [rad]
            phi (_type_): azimuthal angle in a spherical coordinate whose origin is at the halo's center [rad]
        """
        self.zeta = zeta
        self.xi = xi
        self.theta = theta
        self.phi = phi

    def _Norm_factor(self):
        """Normalisation factor

        Returns:
            _type_: _description_
        """
        return (self.zeta * self.xi)**(2.0/3.0)

    @property
    def Cap_A(self):
        """Parameter capital A in Binney 1985

        Returns:
            _type_: _description_
        """
        norm_factor = self._Norm_factor()
        return (np.sin(self.theta)**2 / self.zeta**2 + np.cos(self.theta)**2 / self.xi**2 * (np.sin(self.phi)**2 + np.cos(self.phi)**2 / self.zeta**2)) * norm_factor**2

    @property
    def Cap_B(self):
        """Parameter capital B in Binney 1985

        Returns:
            _type_: _description_
        """
        norm_factor = self._Norm_factor()
        return (np.cos(self.theta) * np.sin(2 * self.phi) * (1 - 1/self.zeta**2) / self.xi**2) * norm_factor**2

    @property
    def Cap_C(self):
        """Parameter capital C in Binney 1985

        Returns:
            _type_: _description_
        """
        norm_factor = self._Norm_factor()
        return ((np.cos(self.phi)**2 + np.sin(self.phi)**2 / self.zeta**2) / self.xi**2) * norm_factor**2

    @property
    def Lcase_f(self):
        """Parameter lowercase f in Binney 1985

        Returns:
            _type_: _description_
        """
        norm_factor = self._Norm_factor()
        return ((np.cos(self.theta)**2 / self.xi**2) + np.sin(self.theta)**2 * (np.cos(self.phi)**2 + (np.sin(self.phi)**2 / self.zeta**2))) * norm_factor

    def AxisRatio(self):
        """Calculate apparent ratio of the projected ellipse from 3d axis ratios and angles of projection

        Returns:
            _type_: axis ratio Q
        """
        A = self.Cap_A
        B = self.Cap_B
        C = self.Cap_C

        numerator = A + C - np.sqrt((A - C)**2 + B**2)
        denominator = A + C + np.sqrt((A - C)**2 + B**2)
        axis_ratio = np.sqrt(numerator / denominator)

        return axis_ratio

    def Ellipticity(self):
        """Ellipticity of the projected ellipse

        Returns:
            _type_: _description_
        """
        axis_ratio = self.AxisRatio()
        ellipticity = (1 - axis_ratio) / (1 + axis_ratio)
        return ellipticity

    def PA_phi(self):
        """compute the position angle of the apparent major axis of the ellipse and the intrinsic major axis of the 3d ellipsoid

        Returns:
            _type_: position angle [rad]
        """
        cap_A = self.Cap_A
        cap_B = self.Cap_B
        cap_C = self.Cap_C

        if cap_A == cap_C:
            return 0.
        else:
            doub_psi = np.arctan(cap_B / (cap_A - cap_C))
            doub_psi = np.where(doub_psi < 0, doub_psi + np.pi, doub_psi)
            return doub_psi / 2

    def RadialProfile(self, R, profile, R_align = 'major', interpolate = False):
        """Compute the radial profile of a projected 3d triaxial halo

        Args:
            R (_type_): radius variable 
            profile (_type_): a "profile" class instance, contains the integrand for projection alopng the LoS
            R_align (str, optional): which principle axis of the projected ellipse the radius varaible is aligned to. Defaults to 'major'. Supported R_align: ['major', 'minor', 'average']. The "average" radial profile is the azimuthal average.
            interpolate (bool, optional): whether to interpolate the radial profile. Defaults to False.

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        cap_A = self.Cap_A
        cap_B = self.Cap_B
        cap_C = self.Cap_C
        lcase_f = self.Lcase_f
        eigenvalue_minus = (cap_A + cap_C - np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2
        eigenvalue_plus = (cap_A + cap_C + np.sqrt((cap_A - cap_C)**2 + cap_B**2)) / 2

        if R_align == 'major':

            as_sq_diag = eigenvalue_minus / lcase_f * R**2
            Sigma_rad = np.zeros_like(as_sq_diag)

            for i, as_sq in enumerate(as_sq_diag):
                Sigma_rad[i] = quad(profile.Project_integrand, 0, np.inf, args = (as_sq))[0] * 2 / np.sqrt(lcase_f)

            if interpolate:
                return interp1d(x=R, y=Sigma_rad, fill_value='extrapolate')
            else:
                return Sigma_rad, R
    
        elif R_align == 'minor':

            as_sq_diag = eigenvalue_plus / lcase_f * R**2
            Sigma_rad = np.zeros_like(as_sq_diag)

            for i, as_sq in enumerate(as_sq_diag):
                Sigma_rad[i] = quad(profile.Project_integrand, 0, np.inf, args = (as_sq))[0] * 2 / np.sqrt(lcase_f)

            if interpolate:
                return interp1d(x=R, y=Sigma_rad, fill_value='extrapolate')
            else:
                return Sigma_rad, R
        
        elif R_align == 'average':

            axis_ratio = self.AxisRatio()
            as_sq_diag = eigenvalue_minus / lcase_f * R**2
            Sigma_rad = np.zeros_like(as_sq_diag)

            for i, as_sq in enumerate(as_sq_diag):
                Sigma_rad[i] = quad(profile.Project_integrand, 0, np.inf, args = (as_sq))[0] * 2 / np.sqrt(lcase_f)

            if interpolate:
                return interp1d(x=R * np.sqrt(axis_ratio), y=Sigma_rad, fill_value='extrapolate')
            else:
                return Sigma_rad, np.sqrt(axis_ratio) * R
        
        else:
            raise ValueError("Supported R_align: ['major', 'average', 'minor'].")


    def Project_2d(self, x_grid, y_grid, profile):
        """calculate the projected light/mass profile on a plane

        Args:
            x_grid (_type_): a grid of x coordinates
            y_grid (_type_): a gird of y coordinates
            profile (_type_): a "profile" class instance, contains the integrand for projection along the LoS

        Returns:
            _type_: projected light/mass profile on a plane
        """
        cap_A = self.Cap_A
        cap_B = self.Cap_B
        cap_C = self.Cap_C
        lcase_f = self.Lcase_f

        as_sq = (cap_A * x_grid**2 + cap_B * x_grid * y_grid + cap_C * y_grid**2) / lcase_f
        Sigma = np.zeros_like(as_sq)

        for i in range(as_sq.shape[0]):
            for j in range(as_sq.shape[1]):
                Sigma[i, j] = quad(profile.Project_integrand, 0, np.inf, args= (as_sq[i, j]))[0] * 2 / np.sqrt(lcase_f)

        return Sigma

    