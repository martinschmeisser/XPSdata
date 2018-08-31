# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 19:03:34 2017

@author: jdo
"""

import numpy as np
import matplotlib.pyplot as plt
from XPS.shirley import remove_shirley_background


def shirley_peak_area(region, e_b_min, e_b_max, doplot=False, name=None):
    """integrate area under peak after shirley background subtraction
    will integrate over whole region, no peak fitting

    Parameters
    ----------
    filename : string
        location of the data file, expects an SpecsLab xml file.
    E_b_min : float
        region boundary
    E_b_max : float
        region boundary
    casa_signal:
        allows to plot a reference signal (default None)
    doplot : bool
        show the signal and background plots (default False)
    """
    # determine index of left and right boundary
    index1 = np.where(abs(region.x_be-(e_b_min)) < 1e-10)[0][0]
    index2 = np.where(abs(region.x_be-(e_b_max)) < 1e-10)[0][0]

    energy = region.x_be                                # binding energy scale
    signal = region.y_avg_cps_mcd                       # signal
    transm = region.info.get('transmission')[0]         # transmission function
#    try:
#        transm = region.info.get('transmission')[0]         # transmission function
#    except:
#        transm = np.ones_like(energy)

    # NOTE
    # original specs signal and casa's export are not scaled for transmission
    # scale s by the transmission function
    signal = signal/transm

    # keep only the peak
    signal = signal[index1:index2]
    energy = energy[index1:index2]
    transm = transm[index1:index2]

    # signal with background subtracted
    _, reduced_signal, background = remove_shirley_background(signal, energy, 30, 1e-6)

    if doplot:
        plt.figure()
        plt.plot(-energy, signal, label='specs signal')  # original signal
        plt.plot(-energy, reduced_signal, label='after shirley subtr.')
        plt.plot(-energy, background, label='shirley bg')  # background

        plt.legend()
        plt.gca().invert_xaxis()
        plt.xlabel('binding energy (eV)')
        plt.ylabel('intensity (cps)')
        fig = plt.gcf()
        fig.savefig(name+'.png', dpi=400)
    # fi

    area = np.trapz(reduced_signal, energy)
    return area, reduced_signal
