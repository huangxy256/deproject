from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

from deproject.Util.orientation import Isotropic_inclination

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
        return np.sqrt((self.qintr * np.sin(inc))**2 + np.cos(inc)**2) # eq.(35) of Cappellari (2020, MNRAS)
    
    @staticmethod
    def _Qobs(qintr, inc):
        """Compute the projected axis ratio

        Args:
            qintr (arr of float): intrinsic axis ratio of axisymmetric ellipsoid
            inc (array of float): inclination angle [rad]

        Returns:
            array of float: observed axis ratio
        """
        return np.sqrt((qintr * np.sin(inc))**2 + np.cos(inc)**2) # eq.(35) of Cappellari (2020, MNRAS)
        
        
    def Ellipticity(self, inc):
        """Calculate the ellipticity of projected ellipse

        Args:
            inc (float): inclination angle [rad]

        Returns:
            arr of float: ellipticity (1-Q)/(1+Q)
        """
        Qobs = self.Qobs(inc)
        # Qobs = np.where(Qobs > 1, 1 / Qobs, Qobs)
        e = (1 - Qobs) / (1 + Qobs)
        return e
    
    @staticmethod
    def _Ellipticity(qintr, inc):
        """Calculate the ellipticity of projected ellipse

        Args:
            qintr (arr of float): intrinsic axis ratio of axisymmetric ellipsoid
            inc (float): inclination angle [rad]

        Returns:
            arr of float: ellipticity (1-Q)/(1+Q)
        """
        Qobs = Projection_axisym._Qobs(qintr, inc)
        # Qobs = np.where(Qobs > 1, 1 / Qobs, Qobs)
        e = (1 - Qobs) / (1 + Qobs)
        return e
    
    def Eobs_dist_isotropic_inc(self, single_proj=1, plot_scatter=False, plot_2dhist=False, bins_2dhist=10):
        """Model the projected ellipticity for a sample under isotropic inclination angle prior

        Args:
            single_proj (int, optional): number of projections for each single intrinsic axis ratio
            plot_scatter (bool, optional): whether to plot scatter plot of eobs vs inclination. Defaults to False.
            plot_2dhist (bool, optional): whether to plot 2d histogram of eobs vs inclination. Defaults to False.
            bins_2dhist (int, optional): bins of the 2d histogram. Defaults to 10.

        Returns:
            array of float: observed ellipticity
            array of float: isotropic inclination angle [rad]
        """
        qintr_expand = np.repeat(self.qintr, single_proj)
        inc_iso = Isotropic_inclination(len(qintr_expand), 1, deg=0)
        eobs = Projection_axisym._Ellipticity(qintr_expand, inc_iso)

        if plot_scatter:
            plt.figure()
            plt.scatter(eobs, inc_iso, marker = '.')
            plt.xlabel('$e$')
            plt.ylabel('$i$ [rad]')
            plt.show()

        if plot_2dhist:
            plt.figure()
            hist, xedges, yedges = np.histogram2d(eobs, inc_iso, bins=bins_2dhist)
            xcoord, ycoord = np.meshgrid(xedges, yedges)
            plt.pcolormesh(xcoord, ycoord, hist.T)
            plt.colorbar(label = 'counts')
            plt.xlabel('$e$')
            plt.ylabel('$i$ [rad]')
            plt.show()

        return eobs, inc_iso
    

    def Cond_pdf_inc_given_eobs(self, eobs, bins_2dhist=10, plot_1dpdf=False, normalize=False, quiet=False, **kwargs_eobs_inc_2dhist):
        eobs_all, inc_iso = self.Eobs_dist_isotropic_inc(**kwargs_eobs_inc_2dhist)
        hist, xedges, yedges = np.histogram2d(eobs_all, inc_iso, bins=bins_2dhist)
        bin_ind = np.digitize(eobs, xedges) - 1

        if not quiet:
            print('bin index = {}' .format(bin_ind))
            print('{:.4f} < e < {:.4f}' .format(xedges[bin_ind], xedges[bin_ind+1]))

        hist_1d = hist[bin_ind, :]

        if normalize:
            hist_1d = hist_1d / np.sum(hist_1d * np.diff(yedges))

        if plot_1dpdf:
            plt.figure()
            plt.bar(yedges[:-1], hist_1d, align='edge', width=np.diff(yedges))
            plt.ylabel('counts/pdf')
            plt.xlabel('$i$ [rad]')
            plt.title(r'$P(i \vert {:.3f} < e < {:.3f})$' .format(xedges[bin_ind], xedges[bin_ind+1]))
            plt.show()

        return yedges, hist_1d