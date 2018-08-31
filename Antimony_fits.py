# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:29:38 2017

@author: jdo
"""

#from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
from fit_sb_3d import fit_sb_3d
from XPS.get_region import get_region


regions = []
results = []
inputs = []
results_tab = np.zeros((13, 11))

dirn = 'D:/EXPERIMENTS/'

regions.append(get_region(dirn + '2015-08-05_P003/P003_XPS.xml', 2, 3))
regions.append(get_region(dirn + '2015-08-11_P004/XPS_P004.xml', 3, 2))
regions.append(get_region(dirn + '2015-08-24_P005/P005.xml', 2, 2))
regions.append(get_region(dirn + '2015-09-16_P006/XPS_P006.xml', 2, 1))
regions.append(get_region(dirn + '2016-05-02_P007/P007.xml', 3, 1))
regions.append(get_region(dirn + '2016-05-09_P008/P008.xml', 2, 1))
regions.append(get_region(dirn + '2016-06-14_P009/P009.xml', 2, 1))
regions.append(get_region(dirn + '2016-06-22_P010/P010.xml', 2, 1))
regions.append(get_region(dirn + '2016-07-11_P011/P011.xml', 2, 1))
regions.append(get_region(dirn + '2016-09-07_P012/2016-09-07_P012.xml', 1, 3))
regions.append(get_region(dirn + '2016-09-27_P013/P013.xml', 2, 1))
regions.append(get_region(dirn + '2017-02-22_P014/2017-02-22_P014.xml', 3, 1))
regions.append(get_region(dirn + '2017-09-11_P015/2017-09-19_P015_Cs-K-Sb.sle', 0, 0))

e_min = [-539, -539, -539, -539, -539, -539, -539, -542, -539, -541, -539, -540, -540]
e_max = [-520, -520, -520, -520, -520, -520, -520, -522, -520, -520, -520, -520, -520]

names = ['P003', 'P004', 'P005', 'P006', 'P007', 'P008', 'P009', 'P010',
         'P011', 'P012', 'P013', 'P014', 'P015']

# %%

for n, r in enumerate(regions):
    inputs.append((r, e_min[n], e_max[n], names[n], True, names[n]))


if __name__ == '__main__':

    # parallel execution, works from console but not in Spyder
    # p = Pool(5)
    # results = p.map(fit_sb_3d, inputs)

    # sequential execution
    results = map(fit_sb_3d, inputs)

    #          P003  4      5     6     7    8      9    10      11    12    13  14  15
#QEs515 = [0.357, 5.2, 0.175, 4.8, 2.36, 1.28, 4.8,  0.25, 2.91,  0.25,  5.2, 10, 10.5]
#QEs532 = [0.357, 5.2, 0.175, 4.8, 1.58, 0.69, 2.59, 0.25, 1.65,  0.01,  3.8, 7.72, 7.23]

    #results_tab[:, 0] = [0.357, 5.2, 0.175, 4.8, 2.28, 1.28, 4.8, 0.25, 2.91,
    #                     0.25, 5.2, 10, 10.5]
# %%
    results_tab[:,0] = [0.357, 5.2, 0.175, 4.8, 1.58, 0.69, 2.59, 0.25, 1.65,  0.01,  3.8, 7.72, 7.23]
# %%
    with open('Sb-fit-results.txt', 'w') as f:
        for n, r in enumerate(results):

            # results_tab[n,0] = QE
            results_tab[n, 1] = r[1]['Sb3d_32_1_amplitude']
            results_tab[n, 2] = r[1]['Sb3d_32_1_sigma']
            results_tab[n, 3] = r[1]['Sb3d_52_1_amplitude']
            results_tab[n, 4] = r[1]['Sb3d_52_1_sigma']
            results_tab[n, 5] = r[1]['Sb3d_32_2_amplitude']
            results_tab[n, 6] = r[1]['Sb3d_32_2_sigma']
            results_tab[n, 7] = r[1]['Sb3d_52_2_amplitude']
            results_tab[n, 8] = r[1]['Sb3d_52_2_sigma']
            results_tab[n, 9] = r[1]['O1s_amplitude']
            results_tab[n, 10] = r[1]['O1s_sigma']

            f.write('%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' %
                    (names[n], results_tab[n, 0],
                     r[1]['Sb3d_32_1_amplitude'], r[1]['Sb3d_32_1_sigma'],
                     r[1]['Sb3d_52_1_amplitude'], r[1]['Sb3d_52_1_sigma'],
                     r[1]['Sb3d_32_2_amplitude'], r[1]['Sb3d_32_2_sigma'],
                     r[1]['Sb3d_52_2_amplitude'], r[1]['Sb3d_52_2_sigma'],
                     r[1]['O1s_amplitude'], r[1]['O1s_sigma']))

    # %%

    cmap = plt.get_cmap('RdYlGn_r')
    plt.figure(figsize=(6,4))
    plt.gca().set_xscale('log')
    plt.scatter(results_tab[:, 3]/results_tab[:, 7], results_tab[:, 0])
                #,color=cmap(results_tab[:, 9]/results_tab[:, 3]*5))

    for i, n in enumerate(names):
        plt.gca().annotate(n, xy=(results_tab[i, 3]/results_tab[i, 7],  results_tab[i, 0]),
                              xytext=(results_tab[i, 3]/results_tab[i, 7]*0.5,  results_tab[i, 0]+0.2))

    plt.xlabel('Sb 3d reacted / metallic ratio')
    plt.ylabel('QE at 532nm (%)')
    #plt.xlim((-0.1, 1.0))
    plt.ylim((-0.1, 8.5))
#    plt.savefig('QE-vs-Sb-ratio.pdf', dpi=600, bbox_inches='tight')

    plt.figure(figsize=(6,4))
    #plt.gca().set_xscale('log')
    plt.scatter(results_tab[:, 9]/results_tab[:, 3], results_tab[:, 0])
                #,color=cmap(results_tab[:, 3]/results_tab[:, 7]/50))

    for i, n in enumerate(names):
        plt.gca().annotate(n, xy=(results_tab[i, 9]/results_tab[i, 3],  results_tab[i, 0]),
                              xytext=(results_tab[i, 9]/results_tab[i, 3]-0.045,  results_tab[i, 0]+0.2) )
    plt.xlabel('O 1s / Sb 3d peak intensity ratio')
    plt.ylabel('QE at 532nm (%)')
    plt.xlim((-0.1, 1.1))
    plt.ylim((-0.1, 8.5))
#    plt.savefig('QE-vs-O-ratio.pdf', dpi=600, bbox_inches='tight')
