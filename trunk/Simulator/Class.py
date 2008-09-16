"""
GC-MS data simulator
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

import numpy

from pyms.Utils.Error import error, stop
from pyms.Utils.Time import time_str_secs
from pyms.Utils.Math import vector_by_step

class GCMS_simulator(object):

    """
    """

    def __init__(self, rt=None, mz=None):

        """
        """

        if rt == None or mz == None:
            error('Retention time and mz specifications must be given')

        if not is_list(rt) or len(rt) != 3:
            error("'rt' specification list must have exactly three elements")

        if not is_list(mz) or len(mz) != 3:
            error("'mz' specification list must have exactly three elements")

        self.__time_list = self.__set_time_list(rt)
        self.__mass_list = self.__set_mass_list(mz)

        n = len(self.__time_list)
        m = len(self.__mass_list)
        self.__intensity_matrix = numpy.zeros((n,m))

    def  __set_time_list(self, rt_spec):

        """
        """

        rt_start = time_str_secs(rt_spec[0])
        rt_stop = time_str_secs(rt_spec[1])
        rt_step = float(rt_spec[2])

        time_list = vector_by_step(rt_start, rt_stop, rt_step) 

        return time_list

    def  __set_mass_list(self, mz_spec):

        """
        """

        mz_start = float(mz_spec[0])
        mz_stop = float(mz_spec[1])
        mz_step = float(mz_spec[2])

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


