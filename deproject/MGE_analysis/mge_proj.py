import numpy as np 
from deproject.Kinematics.mge_misc import sum_gaussian_components

class MGE_Proj(object):

    def __init__(self, peak_intr, sigma, qintr):
        """Class to manage MGE projection

        Args:
            peak_intr (array of float): peaks of the MGEs describing the intrinsic mass density/luminosity profile
            sigma (array of float): sigmas of the MGEs describing the intrinsic mass density/luminosity profile
            qintr (array of float): intrinsic axis ratios of the MGEs describing the intrinsic mass density/luminosity profile; qintr < 1 for oblates, qintr > 1 for prolates
        """
        self.peak_intr = peak_intr
        self.sigma = sigma 
        self.qintr = qintr

    def qobs(self, inc):
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

    def surf(self, inc):
        """Compute peaks of the MGEs describing the surface density/brightness profile

        Args:
            inc (array of float): inclination angle [rad]

        Returns:
            array of float: peaks of the MGEs describing the surface density/brightness profile
        """
        qobs = self.qobs(inc)
        surf = self.peak_intr * np.sqrt(2 * np.pi) * self.sigma * self.qintr / qobs
        return surf

    def mtot(self):
        """Compute the total mass/luminosity by summation of the intrinsic MGEs

        Returns:
            float: total mass/luminosity
        """
        return (self.qintr * self.peak_intr * (np.sqrt(2 * np.pi) * self.sigma)**3).sum()

    def mtot_proj(self, inc):
        """Compute the total mss/luminosity by summation of the projected MGEs

        Args:
            inc (array of float): inclination angle [rad]

        Returns:
            float: total mass/luminosity
        """
        surf = self.surf(inc)
        qobs = self.qobs(inc)

        return (2 * np.pi * surf * self.sigma**2 * qobs).sum()