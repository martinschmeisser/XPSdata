# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 20:37:00 2018

@author: jdo
"""

import logging
import numpy as np
from XPS.cached_property import cached_property

_LOG = logging.getLogger(__name__)
# Make the logger follow the logging setup from the caller
_LOG.addHandler(logging.NullHandler())


class RegionBase(object):
    """Class that represents a region

    The class contains attributes for the items listed in the
    'information_names' class variable.

    Some useful ones are:
    * **name**: The name of the region
    * **ID**: The ID it goes by in the SQLite file
    * **region**: Contains information like, DwellTime, AnalysisMethod,
        ScanDelta, ExcitationEnergy etc.
        Note that these are slightly different than in the SpecsLab 2 XML
        format

    All auxiliary information is also available from the 'info'
    attribute.

    """
    # pylint: disable=too-many-instance-attributes
    # 17 is reasonable in this case.

    name = ''
    region = object
    mcd_head = -1
    mcd_tail = -1
    mcd_num = -1
    detectors = object

    def __repr__(self):
        """Returns class representation"""
        return '<{}(name=\'{}\')>'.format(self.__class__.__name__, self.name)

    @cached_property
    def x(self):  # pylint: disable=invalid-name
        """Returns the kinetic energy x-values as a Numpy array"""
        # Calculate the x-values
        start = self.region['kinetic_energy']
        end = start + (self.region['values_per_curve'] - 1) *\
              self.region['scan_delta']
        data = np.linspace(start, end, self.region['values_per_curve'])
        _LOG.debug('Creating x values from %i to %i in %i steps',
                   start, end, self.region['values_per_curve'])
        return data

    @cached_property
    def x_be(self):
        """Returns the binding energy x-values as a Numpy array"""
        if self.region['analysis_method'] != 'XPS':
            message = "Analysis_method is {}".format(
                self.region['analysis_method'])
            raise NotXPSException(message)

        # Calculate the x binding energy values
        data = self.x - self.region['excitation_energy']
        _LOG.debug('Creating x_be values from %i to %i in %i steps',
                   data.min(), data.max(), data.size)
        return data

    def y_avg_counts(self):
        """ implement """
        raise NotImplementedError("must be implemented by a real instance")

    @cached_property
    def y_avg_cps(self):
        """Returns the average counts per second from a cycle as a Numpy array"""
        try:
            data = self.y_avg_counts / self.region['dwell_time']
            _LOG.debug('Creating %i y_avg_cps values', data.size)
        except TypeError:
            data = None
        return data

    @cached_property
    def y_avg_counts_mcd(self):
        """Returns the average counts from an MCD scan as a Numpy array"""
        data = self.y_avg_counts

        step = self.region['scan_delta']
        mcd_head = self.mcd_head
        mcd_tail = self.mcd_tail
        mcd_num = self.mcd_num
        e_pass = self.region['pass_energy']
        detectors = self.detectors

        # space for single channel data (_sc)
        y_sc = np.zeros(len(self.x)+mcd_head+mcd_tail)
        x_sum = np.zeros(len(self.x)+mcd_head+mcd_tail+1)
        x_sc = np.zeros(len(self.x)+mcd_head+mcd_tail)

        # extend the region.x with steps for mcd_head and mcd_tail
        x_sum[0:mcd_head] = np.arange(np.min(self.x)-step*mcd_head,
                                      np.min(self.x)-0.5*step, step)
        x_sum[mcd_head:(mcd_head+len(self.x))] = self.x
        x_sum[-(mcd_tail+1):] = np.arange(np.max(self.x)+step,
                                          np.max(self.x)+step*(mcd_tail+1.5),
                                          step)

        y_sum = np.zeros(len(self.x))

        # sum up the channels from a multi-channel detector :
        #  the channels are interlaced in a long 1D array
        #  and have a shift in units of e_pass
        for k in range(0, mcd_num):
            # channels are interlaced, extract every 25th number
            y_sc[:] = data[k::mcd_num]
            # make a local copy of the x array, shifted
            x_sc = x_sum + detectors[k]['shift']*e_pass

            for energy in range(0, len(self.x)):
                # IndeX of the first measured point that is larger in e_kin
                # than the x point we want to estimate
                idx = np.where(x_sc[:] > self.x[energy])[0][0]

                if idx == len(y_sc):
                    # sometimes the last point is out of bounds
                    y_sum[energy] = y_sc[idx-1]
                else:
                    # interpolate y value at x[l] from left and right values
                    y_sum[energy] += y_sc[idx-1] + (y_sc[idx] - y_sc[idx-1])/ \
                        (x_sc[idx]-x_sc[idx-1])*(self.x[energy]-x_sc[idx-1])

        _LOG.debug('Creating y_avg_counts values from  scans')
        return y_sum

    @cached_property
    def y_avg_cps_mcd(self):
        """Returns the average counts per second from an MCD scan as a Numpy array"""
        try:
            data = self.y_avg_counts_mcd / self.region['dwell_time']
            _LOG.debug('Creating %i y_avg_cps values', data.size)
        except TypeError:
            data = None
        return data

class NotXPSException(Exception):
    """Exception for trying to interpret non-XPS data as XPS data"""
    pass
