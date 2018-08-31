# -*- coding: utf-8 -*-
""" Copyright 2007-2016 The HyperSpy developers
# -#
# This file is part of  HyperSpy.
# -#
 # HyperSpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# -#
 # HyperSpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# -#
# You should have received a copy of the GNU General Public License
# along with  HyperSpy.  If not, see <http://www.gnu.org/licenses/>.

###
#original file from
#https://github.com/hyperspy/hyperspy/pull/1013/commits/e17925a9456f97ed8534d3d6cffb8d8eb25ef10a
#
#modified for use outside of hyperspy
###
"""
import numpy as np

def remove_shirley_background(signal, energyscale, max_iter=10, eps=1e-6, debug=False):
    """Remove the inelastic background of photoemission SI by the shirley
    iterative method.

    Parameters
    ----------
    signal : array
        photoemission signal
    energyscale : array
        binding energy of the photoemission spectrum
    max_iter : int
        maximum number of iterations
    eps : float
        convergence limit
    """

    cnt = 0
    left = signal[:3].mean()
    right = signal[-3:].mean()
    background = right * 0.9 * np.ones(signal.shape)
    mean_epsilon = 10*eps
    integral = None
    old_integral = None
    while  (mean_epsilon > eps) and (cnt < max_iter):
        if integral is not None:
            old_integral = integral
        reduced_signal = signal - background
        integral = np.cumsum(reduced_signal[::-1], axis=0)[::-1] * energyscale
        background = (left-right)*integral/integral[0] + right
        if old_integral is not None:
            epsilon = np.abs(integral[0] - old_integral[0])
            mean_epsilon = epsilon.mean()
        cnt += 1
    if debug:
        print("shirley %s iterations \t epsilon: %s" % (cnt, mean_epsilon))

    if np.max(reduced_signal) > np.max(signal):
        print("WARNING: please check shirley reduction results! reduced signal\
              is larger than original (zero peak?)")
        reduced_signal = signal-background

    return epsilon, reduced_signal, background
