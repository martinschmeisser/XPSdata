# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 19:03:34 2017

@author: jdo
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from lmfit.models import PseudoVoigtModel
from lmfit.models import GaussianModel
from lmfit.models import DonaichModel

from XPS.get_region import get_region
from XPS.shirley import remove_shirley_background

plt.style.use('ms_diss')

colors = ['#0072bd', '#edb120', '#7e2f8e','#77ac30', '#4dbeee', '#a2142f', '#7f7f7f','#bcbd22', '17becf', '#d95319']

def fitMo3d(filename, groupNo, regionNo, E_b_min, E_b_max, doplot=False):
    """Fit a Moly 3d signal using elemental and oxide components (two for each peak, 3/2 and 5/2).
    A shirley background is subtracted first.
    The elemental Mo components are fittet to the PseudoVoigt line shape.
    Oxide components are fitted to a Gaussian line shape.

    Parameters
    ----------
    filename : string
        location of the data file, expects an SpecsLab xml file.
    groupNo : int
        index of the group of spectra (from 0)
    regionNo : int
        index of the spectrum that should be fitted (from 0)
    E_b_min : float
        region boundary
    E_b_max : float
        region boundary
    doplot : bool
        wether or not to plot the fit and components (default False)
    """

    #load data from SpecsLab xml file and select region
    region = get_region(filename, groupNo, regionNo)

    # determine index of left and right boundary
    index1 = np.where( abs(region.x_be-(E_b_min)) < 1e-10 )
    index2 = np.where( abs(region.x_be-(E_b_max)) < 1e-10 )

    E = region.x_be            [index1[0][0]:index2[0][0]]    #binding energy scale
    s = region.y_avg_counts_mcd[index1[0][0]:index2[0][0]]    #signal

    #signal with background substacted
    epsilon, sb, B = remove_shirley_background(s,E,3,1e-6)

    ######################
    # set up model components and set constraints for the fit
    ######################

    pre11 = 'Mo3d_32_1_'
    mod11 = PseudoVoigtModel(prefix=pre11)
    pars = mod11.make_params()
    pars.add('delta', value=3, min=2, max=4)
    pars[pre11+ 'amplitude'].set(9000, min=500, max=100000)
    pars[pre11+ 'center'   ].set(-231.2, min=-231.4, max=-231)
    pars[pre11+ 'sigma'    ].set(0.3, min=0.1, max=0.4)
    pars[pre11+ 'fraction' ].set(0.5, min=0.1, max=1)

    mod21 = PseudoVoigtModel(prefix='Mo3d_52_1_')
    pars.update( mod21.make_params())
    pars['Mo3d_52_1_amplitude'].set(expr='3/2*Mo3d_32_1_amplitude') #(12000, min=500, max=100000)
    pars['Mo3d_52_1_center'   ].set(expr='Mo3d_32_1_center+delta') #same distance between 3/2 and 5/2 peaks
    #pars['Mo3d_52_1_sigma'    ].set(expr='1.0*Mo3d_32_1_sigma') #all have same sigma value
    pars['Mo3d_52_1_fraction' ].set(expr='1.0*Mo3d_32_1_fraction') #(0.5, min=0.2, max=1)


    pre12 = 'MoO23d_32_2_'
    mod12 = GaussianModel(prefix=pre12)
    pars.update( mod12.make_params())
    pars[pre12+ 'amplitude'].set(2050, min=500, max=100000)
    pars[pre12+ 'center'   ].set(-234.5, min=-236, max=-234)
    pars[pre12+ 'sigma'    ].set(0.3, min=0.2, max=2)

    pre22 = 'MoO23d_52_2_'
    mod22 = GaussianModel(prefix=pre22)
    pars.update( mod22.make_params())
    pars[pre22+ 'amplitude'].set(expr='3/2*MoO23d_32_2_amplitude') #(1050, min=500, max=100000)
    pars[pre22+ 'center'   ].set(expr='MoO23d_32_2_center+delta') #same distance between 3/2 and 5/2 peaks
    pars[pre22+ 'sigma'    ].set(expr='1.0*MoO23d_32_2_sigma') #all have same sigma value


    pre13 = 'MoO23d_32_3_'
    mod13 = GaussianModel(prefix=pre13)
    pars.update( mod13.make_params())
    pars[pre13+ 'amplitude'].set(2050, min=500, max=100000)
    pars[pre13+ 'center'   ].set(-235.5, min=-237, max=-234)
    pars[pre13+ 'sigma'    ].set(0.3, min=0.2, max=2)

    pre23 = 'MoO23d_52_3_'
    mod23 = GaussianModel(prefix=pre23)
    pars.update( mod23.make_params())
    pars[pre23+ 'amplitude'].set(expr='3/2*MoO23d_32_3_amplitude') #(1050, min=500, max=100000)
    pars[pre23+ 'center'   ].set(expr='MoO23d_32_3_center+delta') #same distance between 3/2 and 5/2 peaks
    pars[pre23+ 'sigma'    ].set(expr='1.0*MoO23d_32_3_sigma') #all have same sigma value

    #composite model is a sum of the components
    mod = mod11+mod21+mod12+mod22+mod13+mod23
    #perform the fit
    out  = mod.fit(sb, pars, x=E)
    print(out.fit_report(min_correl=0.25))

    if (doplot):
        fig = plt.figure(None, figsize=(6, 4))
        #plot individual components
        comps = mod.eval_components(params=out.params, x=E)

        labels = ['Mo','Mo','MoO$_3$','MoO$_3$','MoO$_2$','MoO$_2$','MoO$_2$','MoO$_2$']

        for n, key in enumerate(comps):
            if ((n % 2) == 0):
                plt.plot(-E, comps[key]/10**3, color=colors[n], label=labels[n])
            else:
                plt.plot(-E, comps[key]/10**3, color=colors[n-1])
                #plt.fill_between(-E, 0.0, comps[key]/10**3, facecolor=colors[int(i/2)], alpha=0.5)

        #plot measured signal (dots) and best fit (sum of components)
        plt.plot( -E , sb/10**3 , 'k.', label='measurement')
        #plt.plot( -E , s/10**3 , 'k-')
        #plt.plot( -E , B/10**3 , 'k--')
        plt.plot(-E, out.best_fit/10**3, color=colors[1], label='total fit')
        ax = mpl.pyplot.gca()
        ax.set_xlim(241, -E_b_max)
        plt.xlabel('binding energy (eV)')
        plt.ylabel('counts (arb. u.)')

        print(out.best_values['MoO23d_32_2_center']-out.best_values['Mo3d_32_1_center'])
        print(out.best_values['MoO23d_32_3_center']-out.best_values['Mo3d_32_1_center'])

        ## reference : also plot pure Mo signal
        #load data from SpecsLab xml file and select region
        region = get_region(filename,3,2)

        # determine index of left and right boundary
        index1 = np.where( abs(region.x_be-(-235)) < 1e-10 )
        index2 = np.where( abs(region.x_be-(E_b_max)) < 1e-10 )

        E = region.x_be            [index1[0][0]:index2[0][0]]    #binding energy scale
        s = region.y_avg_counts_mcd[index1[0][0]:index2[0][0]]    #signal

        epsilon, sb, B = remove_shirley_background(s,E,3,1e-6)
        plt.plot( -E , sb/(6*10**3) , 'k--', label='pure Mo')
        plt.legend()

        fig.tight_layout()
        #fig.savefig('D:/Martin/diss_svn/data/Mo-shifts/Mo-oxide-shifts.pdf', bbox_inches='tight', pad_inches=0)
        #fig.savefig('Mo-oxide-shifts.pdf', bbox_inches='tight', pad_inches=0)






def fitpureMo3d(filename, groupNo, regionNo, E_b_min, E_b_max, doplot=False):
    """Fit a Moly 3d signal using elemental and oxide components (two for each peak, 3/2 and 5/2).
    A shirley background is subtracted first.
    The elemental Mo components are fittet to the PseudoVoigt line shape.
    Oxide components are fitted to a Gaussian line shape.

    Parameters
    ----------
    filename : string
        location of the data file, expects an SpecsLab xml file.
    groupNo : int
        index of the group of spectra (from 0)
    regionNo : int
        index of the spectrum that should be fitted (from 0)
    E_b_min : float
        region boundary
    E_b_max : float
        region boundary
    doplot : bool
        wether or not to plot the fit and components (default False)
    """

    #load data from SpecsLab xml file and select region
    region = get_region(filename, groupNo, regionNo)

    #Sb = VAMAS.VAMASExperiment('P006_Sb_surveyf.vms');
    #regions = [data[0][2], data[1][1], data[2][1]]
    print(region.x_be)

    # determine index of left and right boundary
    index1 = np.where( abs(region.x_be-(E_b_min)) < 1e-10 )
    index2 = np.where( abs(region.x_be-(E_b_max)) < 1e-10 )

    E = region.x_be            [index1[0][0]:index2[0][0]]    #binding energy scale
    s = region.y_avg_counts_mcd[index1[0][0]:index2[0][0]]    #signal

    #signal with background substacted
    epsilon, sb, B = remove_shirley_background(s,E,3,1e-6)

    ######################
    # set up model components and set constraints for the fit
    ######################

    pre11 = 'Mo3d_32_1_'
    mod11 = PseudoVoigtModel(prefix=pre11)
    pars = mod11.make_params()
    pars.add('delta', value=3, min=2, max=5)
    pars[pre11+ 'amplitude'].set(9000, min=500, max=100000)
    pars[pre11+ 'center'   ].set(-231, min=-231, max=-230)
    pars[pre11+ 'sigma'    ].set(0.3, min=0.2, max=1)
    pars[pre11+ 'fraction' ].set(0.5, min=0.2, max=1)

    mod21 = PseudoVoigtModel(prefix='Mo3d_52_1_')
    pars.update( mod21.make_params())
    pars['Mo3d_52_1_amplitude'].set(expr='3/2*Mo3d_32_1_amplitude') #(12000, min=500, max=100000)
    pars['Mo3d_52_1_center'   ].set(expr='Mo3d_32_1_center+delta') #same distance between 3/2 and 5/2 peaks
    #pars['Mo3d_52_1_sigma'    ].set(expr='1.0*Mo3d_32_1_sigma') #all have same sigma value
    pars['Mo3d_52_1_fraction' ].set(expr='1.0*Mo3d_32_1_fraction') #(0.5, min=0.2, max=1)

    #composite model is a sum of the components
    mod = mod11+mod21
    #perform the fit
    out  = mod.fit(sb, pars, x=E)

    if (doplot):
        plt.figure()
        #plot individual components
        comps = mod.eval_components(params=out.params, x=E)
        i=0
        for key in comps:
            plt.plot(E, comps[key]/10**3, color=colors[i], label=key)
            plt.fill_between(E, 0.0, comps[key]/10**3, facecolor=colors[i], alpha=0.5)
            i = i+1

        #plot measured signal (dots) and best fit (sum of components)
        plt.plot( E , sb/10**3 , 'k.')
        plt.plot( E , s/10**3 , 'k-')
        plt.plot( E , B/10**3 , 'k--')
        plt.plot(E, out.best_fit/10**3, color=colors[1])
        plt.legend()
        ax = mpl.pyplot.gca()
        ax.set_xlim(E_b_min, E_b_max)
        plt.xlabel('binding energy (eV)')
        plt.ylabel('cps ($10^3$)')
    print(out.fit_report(min_correl=0.25))


def fitAsymMo3d(filename, groupNo, regionNo, E_b_min, E_b_max, doplot=False):
    """Fit a Moly 3d signal using elemental and oxide components (two for each peak, 3/2 and 5/2).
    A shirley background is subtracted first.
    The elemental Mo components are fittet to the PseudoVoigt line shape.
    Oxide components are fitted to a Gaussian line shape.

    Parameters
    ----------
    filename : string
        location of the data file, expects an SpecsLab xml file.
    groupNo : int
        index of the group of spectra (from 0)
    regionNo : int
        index of the spectrum that should be fitted (from 0)
    E_b_min : float
        region boundary
    E_b_max : float
        region boundary
    doplot : bool
        wether or not to plot the fit and components (default False)
    """

    #load data from SpecsLab xml file and select region
    region = get_region(filename, groupNo, regionNo)

    #Sb = VAMAS.VAMASExperiment('P006_Sb_surveyf.vms');
    #regions = [data[0][2], data[1][1], data[2][1]]
    print(region.x_be)

    # determine index of left and right boundary
    index1 = np.where( abs(region.x_be-(E_b_min)) < 1e-10 )
    index2 = np.where( abs(region.x_be-(E_b_max)) < 1e-10 )

    E = region.x_be            [index1[0][0]:index2[0][0]]    #binding energy scale
    s = region.y_avg_counts_mcd[index1[0][0]:index2[0][0]]    #signal

    #signal with background substacted
    epsilon, sb, B = remove_shirley_background(s,E,3,1e-6)

    ######################
    # set up model components and set constraints for the fit
    ######################

    pre11 = 'Mo3d_32_1_'
    mod11 = DonaichModel(prefix=pre11)
    pars = mod11.make_params()
    pars.add('delta', value=3, min=2, max=5)
    pars[pre11+ 'amplitude'].set(9000, min=500, max=100000)
    pars[pre11+ 'center'   ].set(-231, min=-231, max=-230)
    pars[pre11+ 'sigma'    ].set(0.3, min=0.2, max=1)
    pars[pre11+ 'gamma' ].set(0.5, min=0.0, max=1)

    mod21 = DonaichModel(prefix='Mo3d_52_1_')
    pars.update( mod21.make_params())
    pars['Mo3d_52_1_amplitude'].set(9000, min=500, max=100000) #(12000, min=500, max=100000)
    pars['Mo3d_52_1_center'   ].set(expr='Mo3d_32_1_center+delta') #same distance between 3/2 and 5/2 peaks
    pars['Mo3d_52_1_sigma'    ].set(0.3, min=0.2, max=1) #all have same sigma value
    pars['Mo3d_52_1_gamma' ].set(0.5, min=0.0, max=1) #(0.5, min=0.2, max=1)

    #composite model is a sum of the components
    mod = mod11+mod21
    #perform the fit
    out  = mod.fit(sb, pars, x=E)

    if (doplot):
        plt.figure()
        #plot individual components
        comps = mod.eval_components(params=out.params, x=E)
        i=0
        for key in comps:
            plt.plot(E, comps[key]/10**3, color=colors[i], label=key)
            plt.fill_between(E, 0.0, comps[key]/10**3, facecolor=colors[i], alpha=0.5)
            i = i+1

        #plot measured signal (dots) and best fit (sum of components)
        plt.plot( E , sb/10**3 , 'k.')
        plt.plot( E , s/10**3 , 'k-')
        plt.plot( E , B/10**3 , 'k--')
        plt.plot(E, out.best_fit/10**3, color=colors[1])
        plt.legend()
        ax = mpl.pyplot.gca()
        ax.set_xlim(E_b_min, E_b_max)
        plt.xlabel('binding energy (eV)')
        plt.ylabel('cps ($10^3$)')
    print(out.fit_report(min_correl=0.25))