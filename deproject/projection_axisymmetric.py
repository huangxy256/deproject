from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

from deproject.Util.orientation import Isotropic_inclination
from deproject.Util.ellipticity import Axis_ratio2ellipticity, Ellipticity2axis_ratio
from deproject.Util.random_var import Draw_from_hist_with_lo_lim

__all__ = ['Projection_axisym']

class Projection_axisym(object):
    def __init__(self, qintr):
        """class to manage axisymmetric ellipsoid projection

        Args:
            qintr (arr of float): the intrinsic axis ratio of axisymetric ellipsoid. If oblate, qintr < 1; if prolate, qintr > 1
        """
        self.qintr = qintr
        qintr_oblate = np.where(qintr > 1, 1 / qintr, qintr)
        self._qintr_oblate = qintr_oblate

    def Qobs(self, inc):
        """Compute the projected axis ratio

        Args:
            inc (array of float): inclination angle [rad]

        Returns:
            array of float: observed axis ratio, Qobs < 1
        """
        return Projection_axisym._Qobs(self._qintr_oblate, inc)
    
    @staticmethod
    def _Qobs(qintr, inc):
        """Compute the projected axis ratio

        Args:
            qintr (arr of float): intrinsic axis ratio of axisymmetric ellipsoid, qintr < 1
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
        return Projection_axisym._Ellipticity(self._qintr_oblate, inc)
    
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
    
    def Eobs_dist_iso_inc(self, single_proj=1, plot_scatter=False, plot_2dhist=False, bins_2dhist_plot=10):
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
        qintr_expand = np.repeat(self._qintr_oblate, single_proj)
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
            hist, xedges, yedges = np.histogram2d(eobs, inc_iso, bins=bins_2dhist_plot)
            xcoord, ycoord = np.meshgrid(xedges, yedges)
            plt.pcolormesh(xcoord, ycoord, hist.T)
            plt.colorbar(label = 'counts')
            plt.xlabel('$e$')
            plt.ylabel('$i$ [rad]')
            plt.show()

        return eobs, inc_iso
    

    def Cond_pdf_inc_given_eobs(self, eobs, bins_2dhist=10, plot_1dpdf=False, normalize=False, quiet=False, **kwargs_eobs_dist_iso_inc):
        eobs_all, inc_iso = self.Eobs_dist_iso_inc(**kwargs_eobs_dist_iso_inc)
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
    

    @staticmethod
    def Inc_min(Qobs, qintr_min=0.05, deg=False):
        """Compute the minimum inclination angle given a projected axis ratio

        Args:
            Qobs (float): projected axis ratio
            qintr_min (float, optional): minimum intrinsic axis ratio. Defaults to 0.
            deg (bool, optional): return in degree. Defaults to False (rad).

        Returns:
            float: minimum inclination angle given a projected axis ratio 
        """
        rad = np.arcsin(np.sqrt((1 - Qobs**2) / (1 - qintr_min**2)))
        if deg:
            return np.degrees(rad)
        else:
            return rad
        
        
    def Draw_inc_from_ellipticity(self, eobs, num, qintr_min=0.05, plot_hist=False, **kwargs_cond_pdf_inc_given_eobs):
        """Draw inclination from the conditional probability P(i|e) for an input e

        Args:
            eobs (float): input observed ellipticity (single value)
            num (int): number of inclination angles to draw from P(i|e)
            qintr_min (float, optional): minimum intrinsic axis ratio of the deprojected ellipsoid. Defaults to 0.05.
            plot_hist (bool, optional): whether to plot the histogram of the drawn inclination angles. Defaults to False.
            **kwargs_cond_pdf_inc_given_eobs: keyword arguments for method Cond_pdf_inc_given_eobs and Eobs_dist_iso_inc

        Returns:
            arr of float: inclination angle [rad]
        """
        bins, hist_inc = self.Cond_pdf_inc_given_eobs(eobs, **kwargs_cond_pdf_inc_given_eobs)
        Qobs = Axis_ratio2ellipticity(eobs)
        inc_lo_lim = Projection_axisym.Inc_min(Qobs, qintr_min, deg=0)
        inc_draw = Draw_from_hist_with_lo_lim(bins, hist_inc, num, inc_lo_lim)

        if plot_hist:
            plt.figure()
            plt.hist(inc_draw, bins=20)
            plt.xlabel('$i$ [rad]')
            plt.ylabel('counts')
            plt.show()

        return inc_draw

    
    def Recover_inclination(self, eobs_arr, bins_2dhist=10, qintr_min=0.05, **kwargs_eobs_dist_iso_inc):
        """_summary_

        Args:
            eobs_arr (arr of float): Array of observed ellipticity.
            bins_2dhist (int, optional): # of bins of the inc vs e 2d histogram. Defaults to 10.
            qintr_min (float, optional): minimum intrinsic axis ratio of the deprojected ellipsoid. Defaults to 0.05.

        Returns:
            arr of float: recivered inclination angle [rad], size matches with input ellipticity
        """
        eobs_all, inc_iso = self.Eobs_dist_iso_inc(**kwargs_eobs_dist_iso_inc)
        hist, xedges, yedges = np.histogram2d(eobs_all, inc_iso, bins=bins_2dhist)
        bin_ind_all = np.digitize(eobs_arr, xedges) - 1
        Qobs_arr = Ellipticity2axis_ratio(eobs_arr)
        inc_lo_lim = Projection_axisym.Inc_min(Qobs_arr, qintr_min, deg=0)

        inc_draw = np.zeros_like(eobs_arr)

        for i, bin_ind in enumerate(bin_ind_all):
            try:
                hist_1d = hist[bin_ind, :]
                inc_draw[i] = Draw_from_hist_with_lo_lim(yedges, hist_1d, 1, inc_lo_lim[i])[0]
            except IndexError:
                inc_draw[i] = np.nan
                print('Ellipticity {:.4f} out of range!' .format(eobs_arr[i]))

        return inc_draw


    @staticmethod
    def Recover_isotropic_inclination(eobs_arr):
        """Recover incliantion angle for a set of ellipticity assuming the inclination angle is isotropic but has a lower limit set by the ellipticity. Difference with method Recover inclination: does not contain intrinsic axis ratio information

        Args:
            eobs_arr (arr of float): input ellipticity

        Returns:
            _type_: reocevered inclination angle [rad]
        """
        Qobs_arr = Ellipticity2axis_ratio(eobs_arr)
        inc_min = Projection_axisym.Inc_min(Qobs_arr)
        inc_iso = Isotropic_inclination(len(eobs_arr))
        ind_redraw = np.where(inc_iso < inc_min)[0]
        for ind in ind_redraw:
            while inc_iso[ind] < inc_min[ind]:
                inc_iso[ind] = Isotropic_inclination(1)[0]
        return inc_iso
