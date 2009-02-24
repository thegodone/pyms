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
        component matrix C and pure mass spectra S, X = C'S + [noise]

    @return: none
    @rtype: NoneType

    @author: Andrew Isaac
    """

    def __init__(self, rt=None, mz=None, nc=100, sigma=0.01):

        """
        @param rt: Retention Time List
        @type rt: ListType
        @param mz: Mass/Charge List with start, stop and step values as ints
        @type mz: ListType
        @param nc: Total number of components
        @type nc: IntType
        @param sigma: Noise level [0,1] (relative to max intensity=1.0)
        @type sigma: FloatType
        """

        if rt == None or mz == None:
            error('Retention time and mz specifications must be given')

        if not is_list(rt) or len(rt) != 3:
            error("'rtl' specification list must have exactly three elements")

        if not is_list(mz) or len(mz) != 3:
            error("'mzl' specification list must have exactly three elements")

        if not is_str(rt[0]) or not is_str(rt[1]):
            error("start/stop time limits must be strings")

        if not is_int(nc) or nc < 1:
            error("'nc' must be an integer greater than 0")

        if not is_number(sigma) or sigma < 0 or sigma > 1:
            error("'sigma' must be a number in the range [0,1]")

        self.__set_time_list(rt)
        self.__set_mass_list(mz)

        np = len(self.__time_list)
        nm = len(self.__mass_list)

        # generate Components and Spectra, mix and add noise
        self.__C = self.__gen_components(nc, np)

        self.__S = self.__gen_spectra(nc, nm)

        Y = self.__gen_mixture(self.__C, self.__S)

        self.__intensity_matrix = self.__add_noise(Y, sigma)

    def get_intensity_matrix(self):

        """
        @summary: Returns the full intensity matrix

        @return: Intensity matrix
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__intensity_matrix)

    def get_time_list(self):

        """
        @summary: Returns the array of time values

        @return: Array of time values in seconds
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__time_list)

    def get_mass_list(self):

        """
        @summary: Returns the m/z list

        @return: List of m/z
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__mass_list)

    def get_components(self):

        """
        @summary: Returns the pure components' chromatogram

        @return: Matrix of component chromatograms
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__C)

    def get_spectra(self):

        """
        @summary: Returns the pure components' spectra

        @return: Matrix of component spectra
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__S)

    def __gen_components(self, nc, np):

        """
        @summary: Generate the matrix C as num. components by num. chrom. points

        @param nc: Total number of components
        @type nc: IntType
        @param np: Total number of chromatographic time points
        @type np: IntType

        @return: Component matrix
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        # create matrix (as an array) of zeroed floats (doubles)
        C = numpy.zeros((nc, np), 'd')

        # TODO: is this efficient?
        for ii in range(nc):
            ch = chromatogram(np)
            for jj in range(np):
                C[ii,jj] = ch[jj]

        return C

    def __gen_spectra(self, nc, nm):

        """
        @summary: Generate the matrix S as num. components by num. mass channels

        @param nc: Total number of components
        @type nc: IntType
        @param nm: Number of m/z channels
        @type nm: IntType

        @return: Spectra matrix
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        # create matrix (as an array) of zeroed floats (doubles)
        S = numpy.zeros((nc, nm), 'd')

        for ii in range(nc):
            m = mass_spectrum(nm)
            for jj in range(nm):
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
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        # X = C*S';
        X = numpy.dot(C.T,S)
        return X

    def __add_noise(self, X, sigma):

        """
        @summary: Generate random noise

        @param X: Real valued Matrix
        @type X: numpy.ndarray
        @param sigma: Noise level [0,1] (relative to max intensity=1.0)
        @type sigma: FloatType

        @return: Input matrix with added noise
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        # Nm = rand(size(X))*sigma;
        # D = X+Nm;
        Nm = numpy.random.uniform(0.0,1.0,X.shape)*sigma
        D = X+Nm
        return D

    def  __set_time_list(self, rt_spec):

        """
        @summary: Set the array of time values

        @param rt_spec: List of start, stop and step times
        @type rt_spec: ListType

        @return: Array of time values in seconds
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        rt_start = time_str_secs(rt_spec[0])
        rt_step = float(rt_spec[2])
        rt_stop = time_str_secs(rt_spec[1]) + rt_step    # include last point

        self.__time_list = vector_by_step(rt_start, rt_stop, rt_step)

    def  __set_mass_list(self, mz_spec):

        """
        @summary: Set the m/z list

        @param mz_spec: List of start, stop and step m/z
        @type mz_spec: ListType

        @return: List of m/z
        @rtype: numpy.ndarray

        @author: Andrew Isaac
        """

        mz_start = float(mz_spec[0])
        mz_step = float(mz_spec[2])
        mz_stop = float(mz_spec[1]) + mz_step    # include last point

        self.__mass_list = vector_by_step(mz_start, mz_stop, mz_step)
