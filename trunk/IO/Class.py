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
from pyms.Utils.Utils import is_str, is_int, is_array
from pyms.Utils.IO import open_for_writing, close_for_writing 

class IonChromatogram:

    """
    @summary: Models ion chromatogram

    An ion chromatogram can be either m/z channel intensities as a
    function of time, or cumulative intensities over all measured
    m/z. In the latter case the ion chromatogram is total ion
    chromatogram (TIC).

    The attribute '_mass' gives the m/z value of the ion chromatogram.
    If _mass = None, the ion chromatogram is TIC.

    @author: Lewis Lee
    @author: Vladimir Likic
    """

    def __init__(self, ia, time_array, mass=None):

        """
        @param ia: A numpy object representing an array of intensity
            values (i.e. ion chromatogram intensity values)
        @type ia: numpy.ndarray
        @param time_array: A numpy object representing an array of
            retention times
        @type time_array: numpy.ndarray
        @param mass: Mass of ion chromatogram (Null if TIC)
        @type mass: IntType
        """

        if not is_array(ia):
            error("'ia' must be a numpy.")

        if not is_array(time_array):
            error("'time_array' must be a numpy.")

        if mass != None and not is_int(mass):
            error("'mass' must be an integer.")

        if ia.size != time_array.size:
            error("Intensity and time arrays differ in length.")

        self._ia = ia
        self._time_array = time_array
        self._mass = mass
        self._time_step = self.__calc_time_step(time_array) 

    def __len__(self):

        """
        @summary: Returns the length of the IonChromatogram object

        @return: Length of ion chromatogram
        @rtype: IntType
        """

        return self._ia.size

    def __calc_time_step(self, time_array):

        """
        @summary: Calculates the time step

        @param time_array: An array of retention times
        @type time_array: numpy.ndarray
        @return: Time step value 
        @rtype: FloatType
        """

        td_list = []
        for ii in range(time_array.size-1):
            td = time_array[ii+1]-time_array[ii]
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

        if ix < 0 or ix > self._ia.size - 1:
            error("index out of bounds")

        return self._ia[ix]

    def get_intensity_array(self):

        """
        @summary: Returns the entire intensity array

        @return: Intensity array
        @rtype: numpy.ndarray
        """

        return self._ia

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

        if ix < 0 or ix > self._time_array.size - 1:
            error("index out of bounds")

        return self._time_array[ix]

    def get_time_array(self):

        """
        @summary: Returns the entire time array

        @return: Time array
        @rtype: numpy.ndarray
        """

        return self._time_array

    def set_intensity_array(self, ia):

        """
        @summary: Sets the value for the intensity array

        @param ia: An array of intensity values
        @type ia: numpy.ndarray
        @return: none
        @rtype: NoneType
        """

        self._ia = ia

    def get_time_step(self):

        """
        @summary: Returns the time step

        @return: Time step
        @rtype: FloatType
        """

        return self._time_step

    def is_tic(self):

        """
        @summary: Returns True if the ion chromatogram is a total
            ion chromatogram (TIC), or False otherwise

        @return: A boolean value indicating if the ion chromatogram
            is a total ion chromatogram (True) or not (False)
        @rtype: BooleanType
        """

        if self._mass == None:
            return True
        else:
            return False

    def write(self, file_name, minutes=False):

        """
        @summary: Writes the ion chromatogram to the specified file

        @param file_name: A string representing the file name to write
            the ion chromatogram to
        @type file_name: StringType
        @param minutes: A boolean value indicating whether time is in
            minutes (True) or seconds (False)
        @type minutes: BooleanType
        @return: none
        @rtype: NoneType
        """

        if not is_str(file_name):
            error("'file_name' must be a string.")

        fp = open_for_writing(file_name)

        ix = 0
        n = self._ia.size

        while ix < n:
            if minutes:
                fp.write("%8.3f  %10d\n" % (self._time_array[ix]/60.0, \
                        self._ia[ix]))
            else:
                fp.write("%8.3f  %10d\n" % (self._time_array[ix], \
                        self._ia[ix]))
            ix = ix + 1

        close_for_writing(fp)

