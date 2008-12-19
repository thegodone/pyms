"""
GC-MS data simulator
Provides a class to simulate mass spectra of mixtures
"""

 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-8 Vladimir Likic                                    #
 #                                                                           #
 #    This program is free software; you can redistribute it and/or modify   #
 #    it under the terms of the GNU General Public License version 2 as      #
 #    published by the Free Software Foundation.                             #
 #                                                                           #
 #    This program is distributed in the hope that it will be useful,        #
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
 #    GNU General Public License for more details.                           #
 #                                                                           #
 #    You should have received a copy of the GNU General Public License      #
 #    along with this program; if not, write to the Free Software            #
 #    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
 #                                                                           #
 #############################################################################

import copy

import math
import numpy
import random

from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import is_str, is_list, is_number, is_int
from pyms.Utils.Time import time_str_secs
from pyms.Utils.Math import vector_by_step

from pyms.Simulator.Functions import chromatogram, mass_spectrum

class GCMS_simulator(object):

    """
    @summary: Simulate mass spectra of mixtures

        Simulates the data matrix X as a linear mixture of the pure
        component matrix C and pure mass spectra S

    @author: Andrew Isaac
    """

    def __init__(self, rt=None, mz=None, num_components=100, sigma=0.01):

        """
        @summary: GCMS_simulator(rt=None, mz=None, num_components=100,
            sigma=0.01)
        @param rt: Retention Time List
        @type rt: List with start and stop as a string and step as an int
        @param mz: Mass/Charge List with start, stop and step values as ints
        @type mz: List with start, stop and step values as ints
        @param num_components: Total number of components
        @type num_components: Positive Integer
        @param sigma: Noise level (relative to max intensity=1.0)
        @type sigma: Positive Real, 0.0 <= sigma <= 1.0
        """

        if rt == None or mz == None:
            error('Retention time and mz specifications must be given')

        if not is_list(rt) or len(rt) != 3:
            error("'rt' specification list must have exactly three elements")

        if not is_list(mz) or len(mz) != 3:
            error("'mz' specification list must have exactly three elements")

        if not is_str(rt[0]) or not is_str(rt[1]):
            error("start/stop time limits must be strings")

        if not is_int(num_components) or num_components < 1:
            error("num_components must be an integer greater than 0")

        if not is_number(sigma) or sigma < 0 or sigma > 1:
            error("sigma must be a number in the range [0,1]")

        self.__time_list = self.__set_time_list(rt)
        self.__mass_list = self.__set_mass_list(mz)

        num_chrom_points = len(self.__time_list)
        num_channels = len(self.__mass_list)

        # generate Components and Spectra, mix and add noise
        C = self.__gen_components(num_components, num_chrom_points)

        S = self.__gen_spectra(num_components, num_channels)

#        self.__intensity_matrix = self.__gen_mixture(C,S)

        Y = self.__gen_mixture(C,S)

        self.__intensity_matrix = self.__addNoise(Y, sigma)

    def  __set_time_list(self, rt_spec):

        """
        """

        rt_start = time_str_secs(rt_spec[0])
        rt_step = float(rt_spec[2])
        rt_stop = time_str_secs(rt_spec[1]) + rt_step    # include last point

        time_list = vector_by_step(rt_start, rt_stop, rt_step)

        return time_list

    def  __set_mass_list(self, mz_spec):

        """
        """

        mz_start = float(mz_spec[0])
        mz_step = float(mz_spec[2])
        mz_stop = float(mz_spec[1]) + mz_step    # include last point

        mass_list = vector_by_step(mz_start, mz_stop, mz_step)

        return mass_list

    def get_intensity_matrix(self):

        """
        @summary: Returns the full intensity matrix

        @return: Intensity matrix
        @rtype: numpy.ndarray
        """

        return copy.deepcopy(self.__intensity_matrix)

    def get_time_list(self):

        """
        @summary: Returns the array of time values

        @return: Array of time values in seconds
        @rtype: numpy.ndarray
        """

        return copy.deepcopy(self.__time_list)

    def get_mass_list(self):

        """
        @summary: Returns the m/z list

        @return: List of m/z
        @rtype: numpy.ndarray
        """

        return copy.deepcopy(self.__mass_list)

    def __gen_components(self, num_components, num_chrom_points):

        """
        @summary: Generate the matrix C as num_components by num_chrom_points
        @param num_components: Total number of components
        @type num_components: Positive Integer
        @param num_chrom_points: Total number of chromatographic time points
        @type num_chrom_points: Positive Integer
        @return: Component matrix
        @type: Real valued Matrix
        """

        # create matrix (as an array) of zeroed floats (doubles)
        C = numpy.zeros((num_components, num_chrom_points), 'd')

        # TODO: is this efficient?
        for ii in range(num_components):
            ch = chromatogram(num_chrom_points)

            for jj in range(num_chrom_points):
                C[ii,jj] = ch[jj]

        return C

    def __gen_spectra(self, num_components, num_channels):

        """
        @summary: Generate the matrix S as num_components by num_channels
        @param num_components: Total number of components
        @type num_components: Positive Integer
        @param num_channels: Number of m/z channels
        @type num_channels: Positive Integer
        @return: Spectra matrix
        @type: Real valued Matrix
        """
        # create matrix (as an array) of zeroed floats (doubles)
        S = numpy.zeros((num_components, num_channels), 'd')

        for ii in range(num_components):
            m = mass_spectrum(num_channels)
            for jj in range(num_channels):
                S[ii,jj] = m[jj]

        return S

    def __gen_mixture(self, C, S):

        """
        @summary: Calculate the data matrix X as a linear mixture of components
        @param C: Component matrix
        @type C: Real valued Matrix
        @param S: Spectra matrix
        @type S: Real valued Matrix
        @return: Mixture matrix
        @type: Real valued Matrix
        """
        # X = C*S';
        X = numpy.dot(C.T,S)
        return X

    def __addNoise(self, X, sigma):
        """
        @summary: Generate random noise
        @param X: Real valued Matrix
        @type X: Real valued Matrix
        @param sigma: Noise level (relative to max intensity=1.0)
        @type sigma: Positive Real, 0.0 <= sigma <= 1.0
        @return: Input matrix with added noise
        @type: Real valued Matrix
        """
        # Nm = rand(size(X))*sigma;
        # D = X+Nm;
        Nm = numpy.random.uniform(0.0,1.0,X.shape)*sigma
        D = X+Nm
        return D

#%
#% Non-negative matrix factorization (works only in matlab)
#%
#
#%'calculating non-negative matrix factorization'
#%[W,H] = nnmf(X,NP);
#%
#%'done'
