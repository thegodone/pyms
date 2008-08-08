"""
Provides classes for pyms.IO wide use
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

import numpy

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str, is_int, is_array, is_list
from pyms.Utils.IO import open_for_writing, close_for_writing 

class IonChromatogram:

    """
    @summary: Models ion chromatogram

    An ion chromatogram is a set of intensities as a function of retention
    time. This can can be either m/z channel intensities (for example, ion
    chromatograms at m/z=65), or cumulative intensities over all measured
    m/z. In the latter case the ion chromatogram is total ion chromatogram
    (TIC).

    The nature of an IonChromatogram object can be revealed by inspecting
    the value of the attribute '__mass'. This is se to the m/z value of the
    ion chromatogram, or to None for TIC.

    @author: Lewis Lee
    @author: Vladimir Likic
    """

    def __init__(self, ia, time_list, mass=None):

        """
        @param ia: Ion chromatogram intensity values
        @type ia: numpy.ndarray
        @param time_list: A list of ion chromatogram retention times
        @type time_list: ListType
        @param mass: Mass of ion chromatogram (Null if TIC)
        @type mass: IntType
        """

        if not is_array(ia):
            error("'ia' must be a numpy.ndarray object")

        if not is_list(time_list):
            error("'time_list' must be a list")

        if mass != None and not is_int(mass):
            error("'mass' must be an integer")

        if ia.size != len(time_list):
            error("Intensity array and time list differ in length")

        self.__ia = ia
        self.__time_list = time_list
        self.__mass = mass
        self.__time_step = self.__calc_time_step(time_list) 

    def __len__(self):

        """
        @summary: Returns the length of the IonChromatogram object

        @return: Length of ion chromatogram
        @rtype: IntType
        """

        return self.__ia.size

    def __calc_time_step(self, time_list):

        """
        @summary: Calculates the time step

        @param time_list: A list of retention times
        @type time_list: ListType

        @return: Time step value 
        @rtype: FloatType
        """

        td_list = []
        for ii in range(len(time_list)-1):
            td = time_list[ii+1]-time_list[ii]
            td_list.append(td) 

        td_array = numpy.array(td_list)
        time_step = td_array.mean()

        return time_step

    def get_intensity_at_index(self, ix):

        """
        @summary: Returns intensity at given index

        @param ix: An index
        @type ix: IntType

        @return: Intensity value
        @rtype: FloatType
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > self.__ia.size - 1:
            error("index out of bounds")

        return self.__ia[ix]

    def get_intensity_array(self):

        """
        @summary: Returns the entire intensity array

        @return: Intensity array
        @rtype: numpy.ndarray
        """

        return self.__ia

    def get_time_at_index(self, ix):

        """
        @summary: Returns time at given index

        @param ix: An index
        @type ix: IntType

        @return: Time value
        @rtype: FloatType
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > len(self.__time_list) - 1:
            error("index out of bounds")

        return self.__time_list[ix]

    def get_time_list(self):

        """
        @summary: Returns the time list 

        @return: Time list 
        @rtype: ListType
        """

        return self.__time_list

    def set_intensity_array(self, ia):

        """
        @summary: Sets the value for the intensity array

        @param ia: An array of intensity values
        @type ia: numpy.ndarray

        @return: none
        @rtype: NoneType
        """

        self.__ia = ia

    def get_time_step(self):

        """
        @summary: Returns the time step

        @return: Time step
        @rtype: FloatType
        """

        return self.__time_step

    def is_tic(self):

        """
        @summary: Returns True if the ion chromatogram is a total ion
            chromatogram (TIC), or False otherwise

        @return: A boolean value indicating if the ion chromatogram
            is a total ion chromatogram (True) or not (False)
        @rtype: BooleanType
        """

        if self.__mass == None:
            return True
        else:
            return False

    def write(self, file_name, minutes=False):

        """
        @summary: Writes the ion chromatogram to the specified file

        @param file_name: File for writing the ion chromatogram
        @type file_name: StringType
        @param minutes: A boolean value indicating whether to write
            time in minutes
        @type minutes: BooleanType

        @return: none
        @rtype: NoneType
        """

        if not is_str(file_name):
            error("'file_name' must be a string")

        fp = open_for_writing(file_name)

        time_list = self.__time_list

        if minutes:
            for ii in range(len(time_list)):
                time_list[ii] = time_list[ii]/60.0

        for ii in range(len(time_list)):
            fp.write("%8.4f %16.4f\n" % (time_list[ii], self.__ia[ii]))

        close_for_writing(fp)

