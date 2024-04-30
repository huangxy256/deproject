import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def Draw_from_pdf(x_pdf, y_pdf, num, plot_cdf=0, plot_hist=0, bins_hist=10):
    """draw random variable from a probability density function

    Args:
        x_pdf (array of float): x coordinate in the pdf
        y_pdf (array of float): y coordinate of the pdf
        num (int): number of points to draw
        plot_cdf (bool, optional): whether to plot the CDF used to draw the random variable. Defaults to 0.
        plot_hist (bool, optional): whether to plot the histogram. Defaults to 0.
        bins_hist (int, optional): number of bins in the plotted histogram. Defaults to 10.

    Returns:
        array of float: points drawn from the PDF
    """
    cdf = np.cumsum(y_pdf) / np.cumsum(y_pdf)[-1]
    inv_cdf = interp1d(x=cdf, y=x_pdf, fill_value='extrapolate')
    random_num = np.random.uniform(0, 1, num)
    var_rm = inv_cdf(random_num)

    if plot_cdf:
        plt.figure()
        plt.plot(x_pdf, cdf, ls = '-', marker = '.')
        plt.xlabel('x')
        plt.ylabel('CDF')
        plt.show()

    if plot_hist:
        plt.figure()
        plt.hist(var_rm, density=True, bins=bins_hist, histtype='step', lw = 2)
        plt.plot(x_pdf, y_pdf, ls = '-', marker = ' ', lw = 2)
        plt.xlabel('x')
        plt.ylabel('PDF')

    return var_rm


def Draw_from_hist(bins, counts, num):
    """Draw points from a histogram

    Args:
        bins (array of float): bin edges of the histogram
        counts (array of float): counts/density of each bin
        num (int): number of data points to draw

    Returns:
        array of float: random variables drawn from a histogram
    """
    x = np.zeros(num)
    weights = counts / np.sum(counts)
    bin_index = np.random.choice(len(bins)-1, p=weights, size=num)
    for i, bin_ind in enumerate(bin_index):
        bin_start, bin_end = bins[bin_ind], bins[bin_ind+1]
        x[i] = np.random.uniform(bin_start, bin_end)
    return x


def Draw_from_hist_with_lo_lim(bins, counts, num, lower_lim):
    """Draw random variables from a histogram with lower cutoff (reject when lower than the cutoff)

    Args:
        bins (array of float): bin edges of the histogram
        counts (array of float): counts/density of each bin
        num (int): number of data points to draw
        lower_lim (float): lower cutoff

    Returns:
        array of float: random variables drawn from a histogram
    """
    x = Draw_from_hist(bins, counts, num)
    x_ind = np.where(x <= lower_lim)[0]
    for ind in x_ind:
        while x[ind] <= lower_lim:
            x[ind] = Draw_from_hist(bins, counts, 1)[0]
    return x