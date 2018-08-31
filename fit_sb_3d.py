# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 19:03:34 2017

@author: jdo
"""

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import PseudoVoigtModel
from lmfit.models import GaussianModel

from XPS.shirley import remove_shirley_background

COLORS = ['#0072bd', '#d95319', '#edb120', '#7e2f8e', '#77ac30', '#4dbeee',
          '#a2142f', '#7f7f7f', '#bcbd22', '17becf']

# pylint: disable-msg=too-many-statements
# pylint: disable-msg=R0914


def fit_sb_3d(inputs):
    """Fit an Antimony 3d signal using four Sb components (two for each peak,
    3/2 and 5/2). A shirley background is subtracted first.
    The Sb components are fittet to the PseudoVoigt line shape.
    An Oxygen component is fitted to a Gaussian line shape.

    Parameters
    ----------
    region : object
        region class object of the spectrum that should be fitted (from 0)
    e_b_min : float
        region boundary
    e_b_max : float
        region boundary
    doplot : bool
        wether or not to plot the fit and components (default False)
    """

    (region, e_b_min, e_b_max, name, doplot, title) = inputs
    print((region, e_b_min, e_b_max))

    # determine index of left and right boundary
    index1 = np.where(abs(region.x_be-(e_b_min)) < 1e-10)
    index2 = np.where(abs(region.x_be-(e_b_max)) < 1e-10)

    # select energy scale and signal within boundaries
    energy = region.x_be[index1[0][0]:index2[0][0]]
    signal = region.y_avg_counts_mcd[index1[0][0]:index2[0][0]]

    # signal with background substacted
    _, signal_b, background = remove_shirley_background(signal,
                                                        energy, 3, 1e-6)

    ######################
    # set up model components and set constraints for the fit
    ######################

    pre11 = 'Sb3d_32_1_'
    mod11 = PseudoVoigtModel(prefix=pre11)
    pars = mod11.make_params()
    pars.add('delta', value=10, min=9, max=11)          # spin-split value, should be 9.4 eV
    pars.add('delta2', value=-1, min=-2.2, max=-0.9)    # chemical shift between the two Sb species
    pars[pre11 + 'amplitude'].set(9000,    min=0, max=100000)
    pars[pre11 + 'center'   ].set(-535.65, min=-536, max=-535.5)
    pars[pre11 + 'sigma'    ].set(0.3,     min=0.2, max=1)
    pars[pre11 + 'fraction' ].set(0.5,     min=0.2, max=1)
    pre12 = 'Sb3d_32_2_'
    mod12 = PseudoVoigtModel(prefix=pre12)
    pars.update(mod12.make_params())
    pars[pre12 + 'amplitude'].set(2050, min=0, max=100000)
    pars[pre12 + 'center'   ].set(expr='Sb3d_32_1_center+delta2')
    # all have same sigma value
    pars[pre12 + 'sigma'    ].set(expr='1.0*Sb3d_32_1_sigma')
    pars[pre12 + 'fraction' ].set(expr='1.0*Sb3d_32_1_fraction')

    mod21 = PseudoVoigtModel(prefix='Sb3d_52_1_')
    pars.update(mod21.make_params())
    pars['Sb3d_52_1_amplitude'].set(expr='3/2*Sb3d_32_1_amplitude')
    # same distance between 3/2 and 5/2 peaks
    pars['Sb3d_52_1_center'   ].set(expr='Sb3d_32_1_center+delta')
    # all have same sigma value
    # pars['Sb3d_52_1_sigma'    ].set(expr='1.0*Sb3d_32_1_sigma')
    pars['Sb3d_52_1_fraction' ].set(expr='1.0*Sb3d_32_1_fraction')
    mod22 = PseudoVoigtModel(prefix='Sb3d_52_2_')
    pars.update(mod22.make_params())
    pars['Sb3d_52_2_amplitude'].set(expr='3/2*Sb3d_32_2_amplitude')
    # same distance between 3/2 and 5/2 peaks
    pars['Sb3d_52_2_center'].set(expr='Sb3d_32_2_center+delta')
    # all have same sigma value
    # pars['Sb3d_52_2_sigma'    ].set(expr='1.0*Sb3d_32_1_sigma')
    pars['Sb3d_52_2_fraction'].set(expr='1.0*Sb3d_32_1_fraction')

    mod3 = GaussianModel(prefix='O1s_')
    pars.update(mod3.make_params())
    pars['O1s_amplitude'].set(2050,   min=0, max=500000)
    pars['O1s_center'   ].set(-531.9, min=-532.5, max=-531)
    pars['O1s_sigma'    ].set(1.3,    min=0.2, max=2)

    # composite model is a sum of the components
    mod = mod11 + mod12 + mod21 + mod22 + mod3
    # perform the fit
    out = mod.fit(signal_b, pars, x=energy)

    if doplot:
        fig = plt.figure(figsize=(6,4))
        # plot individual components
        comps = mod.eval_components(params=out.params, x=energy)
        labels = ('Sb 3d 3/2 1', 'Sb 3d 3/2 2', 'Sb 3d 5/2 1', 'Sb 3d 5/2 2', 'O 1s')
        for i, key in enumerate(comps):
            plt.plot(-energy, comps[key]/10**3, color=COLORS[i], label=labels[i])
            plt.fill_between(-energy, 0.0, comps[key]/10**3,
                             facecolor=COLORS[i], alpha=0.5)

        # plot measured signal (dots) and best fit (sum of components)
        plt.plot(-energy, signal_b/10**3, 'k.')
        plt.plot(-energy, signal/10**3, 'k-')
        plt.plot(-energy, background/10**3, 'k--')
        plt.plot(-energy, out.best_fit/10**3, color=COLORS[1])
        plt.legend()
        ax = plt.gca()
        ax.set_xlim(540, 520)
        ax.xaxis.set_ticks(np.arange(519, 540, 5.0))
        plt.xlabel('binding energy (eV)')
        plt.ylabel('counts ($10^3$/s)')

        if title is not None:
            plt.title(title)
            plt.savefig('Sb3d-fits/Sb3d-fit-%s.pdf' % title, dpi=600, bbox_inches='tight')
            #plt.close(fig)
    # fi

#       ### figures below are for debugging ###
#    plt.figure(33)
#    plt.plot(energy, signal_b/np.max(signal_b))
#    plt.figure(34)
#    plt.plot(energy, out.best_fit/np.max(out.best_fit))
#    plt.figure(35)
#    plt.plot(energy, signal_b/np.max(signal_b) - out.best_fit/np.max(out.best_fit))
#
#    comps = mod.eval_components(params=out.params, x=energy)
#    plt.figure(36)
#    a = comps['Sb3d_32_1_']
#    b = comps['Sb3d_52_1_']
#    c = comps['Sb3d_32_2_']
#    d = comps['Sb3d_52_2_']
#    maxpeak = np.max((np.amax(a), np.amax(b), np.amax(c), np.amax(d)))
#    plt.plot(energy, (a+b) / maxpeak)
#
#    plt.figure(37)
#    plt.plot(energy, (c+d) / maxpeak)

    # print(out.fit_report(min_correl=0.25))
    return (name, out.best_values, np.sum(np.abs(signal_b-out.best_fit)))
