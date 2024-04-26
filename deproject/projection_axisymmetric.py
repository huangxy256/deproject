from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np

__all__ = ['Projection_axisym']

class Projection_axisym(object):
    def __init__(self, qintr):
        """class to manage axisymmetric ellipsoid projection

        Args:
            qintr (arr of float): the intrinsic axis ratio of axisymetric ellipsoid. If oblate, qintr < 1; if prolate, qintr > 1
        """
        self.qintr = qintr

    def Qobs(self, inc):
        """Compute the projected axis ratio

        Args:
            inc (array of float): inclination angle [rad]

        Returns:
            array of float: observed axis ratio
        """
        if self.qintr <= 1:
            return np.sqrt((self.qintr * np.sin(inc))**2 + np.cos(inc)**2) # eq.(35) of Cappellari (2020, MNRAS)
        else:
            return self.qintr / np.sqrt(np.sin(inc)**2 + self.qintr**2 * np.cos(inc)**2) # eq.(36) of Cappellari (2020, MNRAS)
        
    def Ellipticity(self, inc):
        """Calculate the ellipticity of projected ellipse

        Args:
            inc (float): inclination angle [rad]

        Returns:
            arr of float: ellipticity (1-Q)/(1+Q)
        """
        Qobs = self.Qobs(inc)
        Qobs = np.where(Qobs > 1, 1 / Qobs, Qobs)
        e = (1 - Qobs) / (1 + Qobs)
        return e