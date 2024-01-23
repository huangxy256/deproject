import matplotlib.pyplot as plt 
import numpy as np

def sum_gaussian_components(x, peaks, sigmas):
    """sum up the total counts of the Gaussian components

    Args:
        x (_type_): _description_
        peaks (_type_): _description_
        sigmas (_type_): _description_

    Returns:
        float: total counts of all the Gaussians
    """
    total = np.zeros_like(x)

    for i, r in enumerate(x):
        total[i] = np.sum(peaks * np.exp(-1 * r**2 / 2 / sigmas**2))

    return total

def plot_mge(amplitude, sigma, rad, radial_profile=None, plot_radial_profile=0):

    fig = plt.figure(figsize = (7, 7))
    gs = fig.add_gridspec(2, hspace = 0.03)
    ax = gs.subplots(sharex = True)

    gs_sum = sum_gaussian_components(rad, amplitude, sigma)
    if plot_radial_profile:
        ax[0].plot(rad, radial_profile, marker = 'o', label = 'radial profile')
        ax[1].plot(rad, (radial_profile - gs_sum) / radial_profile * 100., c = ax[0].get_lines()[1].get_color(), marker = ' ')
    ax[0].plot(rad, gs_sum, ls = '-', marker = ' ', label = 'Gaussian sum')

    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].set_ylabel('counts')

    ax[0].legend()

    ax[1].set_ylim([-15, 15])
    ax[1].axhline(0., ls = '--', c = ax[0].get_lines()[0].get_color())
    ax[1].set_ylabel('residuals [%]')
    ax[1].set_xlabel('R [arcsec]')
    