# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:29:38 2017

@author: jdo
"""

import numpy as np
import matplotlib.pyplot as plt
from XPS.get_region import get_region
from XPS.shirley_peak_area import shirley_peak_area
# ## pylint: disable=locally-disabled, multiple-statements, fixme, line-too-long
# pylint: disable=line-too-long


def CsKSb_stoichiometry(region):
    """ calculate the stoichiometry from a CsKSb survey sprectrum
        using hard-coded IMFPs and cross sections """

    # Choose inelastic mean free path and cross sections
    # depending on excitation energy Al Ka or Mg Ka.
    # Numbers from Julius' xlsx.
    # if region.source_info['excitations'][0]['energy'] == 1486.61:
    if region.region['excitation_energy'] == 1486.61:  # AlKa source
        # print('Al Ka excitation\n')
        # E_hv = 1486.61
        k2p_imfp = 41.157
        cs3d_imfp = 28.994
        sb3d_imfp = 34.684

        k2p_cs = 3.97
        cs3d_cs = 40.2
        sb3d_cs = 27.7
    elif region.region['excitation_energy'] == 1253.6:  # MgKa source
        # print('Mg Ka excitation\n')
        # E_hv = 1253.6
        k2p_imfp = 34.701
        cs3d_imfp = 21.998
        sb3d_imfp = 27.979

        k2p_cs = 4.04
        cs3d_cs = 38.73
        sb3d_cs = 27.26
    # fi

    cs3d_area, _ = shirley_peak_area(region, -751, -720, False, 'Cs3d')
    cs3d_conc = cs3d_area/(cs3d_imfp*cs3d_cs)
    # print('Cs 3d %11.2f \t  %11.2f' % (cs3d_area, cs3d_conc))

    k2p_area, _ = shirley_peak_area(region, -309, -288, False, 'K2p')
    k2p_conc = k2p_area/(k2p_imfp*k2p_cs)
    # print('K  2p %11.2f \t  %11.2f' % (k2p_area, k2p_conc))

    sb3d_area, _ = shirley_peak_area(region, -548, -521, False, 'Sb3d')
    sb3d_conc = sb3d_area/(sb3d_imfp*sb3d_cs)
    # print('Sb 3d %11.2f \t  %11.2f\n' % (sb3d_area, sb3d_conc))

    sum_ = k2p_area + sb3d_area + cs3d_area
    cs3d_area = cs3d_area / sum_
    k2p_area = k2p_area / sum_
    sb3d_area = sb3d_area / sum_
    print('Cs K Sb : %.3f \t  %.3f \t %.3f \t relative peak areas' %
          (cs3d_area, k2p_area, sb3d_area))

    sum_ = k2p_conc + sb3d_conc + cs3d_conc
    cs3d_conc = cs3d_conc / sum_
    k2p_conc = k2p_conc / sum_
    sb3d_conc = sb3d_conc / sum_
    print('Cs K Sb : %.3f \t  %.3f \t %.3f \t relative peak intensities (scaled by IMFP)' %
          (cs3d_conc, k2p_conc, sb3d_conc))

#    k2p_conc = k2p_conc / sb3d_conc
#    cs3d_conc = cs3d_conc / sb3d_conc
#    sb3d_conc = sb3d_conc / sb3d_conc
#    print('Cs K Sb : %.3f \t  %.3f \t %.3f \t stoichiometry' %
#          (cs3d_conc, k2p_conc, sb3d_conc))

    return ((cs3d_area, k2p_area, sb3d_area), (cs3d_conc, k2p_conc, sb3d_conc),
            (cs3d_conc / sb3d_conc, k2p_conc / sb3d_conc, sb3d_conc / sb3d_conc))
# end CsKSb_stoichiometry


if __name__ == '__main__':

    filenames = ['D:/EXPERIMENTS/2015-08-05_P003/P003_XPS.xml',
                 'D:/EXPERIMENTS/2015-08-11_P004/XPS_P004.xml',
                 'D:/EXPERIMENTS/2015-08-24_P005/P005.xml',
                 'D:/EXPERIMENTS/2015-09-16_P006/XPS_P006.xml',
                 'D:/EXPERIMENTS/2016-05-02_P007/P007.xml',
                 'D:/EXPERIMENTS/2016-05-09_P008/P008.xml',
                 'D:/EXPERIMENTS/2016-06-14_P009/P009.xml',
                 'D:/EXPERIMENTS/2016-06-22_P010/P010.xml',
                 'D:/EXPERIMENTS/2016-07-11_P011/P011.xml',
                 'D:/EXPERIMENTS/2016-09-07_P012/2016-09-07_P012.xml',
                 'D:/EXPERIMENTS/2016-09-27_P013/P013.xml',
                 'D:/EXPERIMENTS/2017-02-22_P014/2017-02-22_P014.xml',
                 'D:/EXPERIMENTS/2017-09-11_P015/2017-09-19_P015_Cs-K-Sb.sle']

    groupNos = [2,  # P003
                3,  # P004
                2,  # P005
                2,  # P006
                3,  # P007
                2,  # P008
                2,  # P009
                2,  # P010
                2,  # P011
                1,  # P012
                2,  # P013
                3,  # P014
                1]  # P015

    regionNos = [0,  # P003
                 1,  # P004
                 1,  # P005
                 0,  # P006
                 0,  # P007
                 0,  # P008
                 0,  # P009
                 0,  # P010
                 4,  # P011
                 2,  # P012
                 4,  # P013
                 0,  # P014
                 0]  # P015

    ## %%
    stoich = np.zeros((len(filenames),  3))

    for n, filename in enumerate(filenames):
        reg = get_region(filename, groupNos[n], regionNos[n])
        areas, ints, stoich[n, :] = CsKSb_stoichiometry(reg)

        print(areas)
        print(ints)
        print(stoich)

    print('Cs avg.\t', np.mean(stoich[:, 0]))
    print('K avg.\t', np.mean(stoich[:, 1]))
    print('Cs+K avg.\t', np.mean(stoich[:, 0]+stoich[:, 1]))

    # %%
    #                  P003     4     5     6     7    8     9    10    11    12      13   14   15
    QEs515 = np.array([0.357, 5.2, 0.175, 4.8, 2.36, 1.28, 4.8,  0.25, 2.91,  0.25,  5.2, 10,   10.5])
    QEs532 = np.array([0.357, 5.2, 0.175, 4.8, 1.58, 0.69, 2.59, 0.25, 1.65,  0.01,  3.8, 7.72, 7.23])

    names = ['P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010',
             'P011', 'P012', 'P013', 'P014', 'P015']

    results_tab = np.loadtxt('Sb-fit-results.txt', usecols=range(1, 12))
    cmap = plt.get_cmap('RdYlGn_r')
    colors = cmap(results_tab[:, 9]/results_tab[:, 3]*5)


    # ###############################
    # prepare the plots for 532 nm
    # ###############################

    #names = ['P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010',  'P011', 'P012', 'P013', 'P014', 'P015']
    index = [  0,      1,      2,      3,      4,       5,     6,      7,       8,      9,      10,     11,      12]
    #index =  [          1,              3,      4,       5,     6,               8,                      11,      12]

    fig, ax = plt.subplots(figsize=(6.5,3.5))
    ax.scatter(stoich[index, 0], QEs532[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 0], QEs532[i]),
                              xytext=(stoich[i, 0]-0.1, QEs532[i]+0.2) )
    ax.set_xlabel('Cs stoichiometric content')
    ax.set_ylabel('QE (%) at 532 nm')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=0.2))
    # fake up the array of the scalar mappable...
    # sm._A = []
    sm.set_array([])
    cb = fig.colorbar(sm, ticks=(0, 0.2), extend='max')
    cb.set_label('rel. oxygen signal')
    fig.tight_layout()
    plt.xlim((0.0, 3.0))
    plt.ylim((-0.1, 8.5))
    fig.savefig('QE-vs-Cs-content-532nm-all.pdf', bbox_inches='tight', pad_inches=0)


    fig = plt.figure(figsize=(5.5,3.5))
    plt.scatter(stoich[index, 1], QEs532[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 1], QEs532[i]),
                              xytext=(stoich[i, 1]-0.1, QEs532[i]+0.2) )
    plt.xlabel('K stoichiometric content')
    plt.ylabel('QE (%) at 532 nm')
    plt.xlim((0.0, 3.0))
    plt.ylim((-0.1, 8.5))
    fig.tight_layout()
    fig.savefig('QE-vs-K-content-532nm-all.pdf', bbox_inches='tight', pad_inches=0)

    fig = plt.figure(figsize=(6,3.5))
    plt.scatter(stoich[index, 0] + stoich[index, 1], QEs532[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 0] + stoich[i, 1], QEs532[i]),
                              xytext=(stoich[i, 0] + stoich[i, 1]-0.1, QEs532[i]+0.2) )
    plt.xlabel('Cs+K stoichiometric content')
    plt.ylabel('QE (%) at 532 nm')
    plt.xlim((-0.5, 6.0))
    plt.ylim((-0.1, 8.5))
    fig.tight_layout()
    fig.savefig('QE-vs-Alkali-content-532nm-all.pdf', bbox_inches='tight', pad_inches=0)


    fig = plt.figure(figsize=(6,4))
    plt.scatter(stoich[index, 0],stoich[index, 1], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 0], stoich[i, 1]),
                              xytext=(stoich[i, 0]-0.1, stoich[i, 1]+0.2) )
    plt.xlabel('Cs stoichiometric content')
    plt.ylabel('K stoichiometric content')
    plt.xlim(0,3)
    plt.ylim(0,3)
    plt.plot([0, 3],[3, 0], 'k--')
    fig.tight_layout()
    fig.savefig('Cs-vs-K-content-all.pdf', bbox_inches='tight', pad_inches=0)


    # ###############################
    # prepare the plots for 515 nm
    # ###############################

    #names = ['P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010',  'P011', 'P012', 'P013', 'P014', 'P015']
    index = [                                   4,       5,     6,      7,       8,      9,      10,     11,      12]
    #index =  [                                  4,       5,     6,               8,                      11,      12]

    fig, ax = plt.subplots(figsize=(6.5,3.5))
    ax.scatter(stoich[index, 0], QEs515[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 0], QEs515[i]),
                              xytext=(stoich[i, 0]-0.1, QEs515[i]+0.2) )
    ax.set_xlabel('Cs stoichiometric content')
    ax.set_ylabel('QE (%) at 515 nm')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=0.2))
    # fake up the array of the scalar mappable...
    # sm._A = []
    sm.set_array([])
    cb = fig.colorbar(sm, ticks=(0, 0.2), extend='max')
    cb.set_label('rel. oxygen signal')
    fig.tight_layout()
    plt.xlim((0.0, 3.0))
    plt.ylim((-0.1, 11))
    fig.savefig('QE-vs-Cs-content-515nm-all.pdf', bbox_inches='tight', pad_inches=0)


    fig = plt.figure(figsize=(5.5,3.5))
    plt.scatter(stoich[index, 1], QEs515[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 1], QEs515[i]),
                              xytext=(stoich[i, 1]-0.1, QEs515[i]+0.2) )
    plt.xlabel('K stoichiometric content')
    plt.ylabel('QE (%) at 515 nm')
    plt.xlim((-0.5, 3.0))
    plt.ylim((-0.1, 11))
    fig.tight_layout()
    fig.savefig('QE-vs-K-content-515nm-all.pdf', bbox_inches='tight', pad_inches=0)


    fig = plt.figure(figsize=(6,4))
    plt.scatter(stoich[index, 0] + stoich[index, 1], QEs515[index], color=colors[index], edgecolors='k')
    for i in index:
        plt.gca().annotate(names[i], xy=(stoich[i, 0] + stoich[i, 1], QEs515[i]),
                              xytext=(stoich[i, 0] + stoich[i, 1]-0.1, QEs515[i]+0.2) )
    plt.xlabel('Cs+K stoichiometric content')
    plt.ylabel('QE (%) at 515 nm')
    plt.xlim((-0.5, 6.0))
    plt.ylim((-0.1, 11))
    fig.tight_layout()
    fig.savefig('QE-vs-Alkali-content-515nm-all.pdf', bbox_inches='tight', pad_inches=0)


# %%
#    index =  [          1,              3,      4,       5,     6,               8,                      11,      12]
#    fig = plt.figure(figsize=(6,4))
#    plt.scatter(stoich[index, 0], stoich[index, 1])
#    plt.xlim(0,3)
#    plt.ylim(0,3)
#    plt.plot([0, 3],[3, 0])

#end if main