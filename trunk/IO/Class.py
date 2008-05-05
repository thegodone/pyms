"""Class.py
 Module Class in metab.IO
 Contains classes used in reading raw data.
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

from metab.Utils.Error import error
from metab.Utils.Utils import * 
from metab.Utils.IO import open_for_writing, close_for_writing 

class IonChromatogram:

    """class IonChromatogram

    Models IonChromatogram.

    __init__(self, ia, time_array, mass=None)

    @param ia A numpy object representing an
        array of intensity values (i.e. ion chromatogram).
    @param time_array A numpy object representing an array of the
        retention times.
    @param mass An integer value representing the mass of the
        ion chromatogram.
    """

    def __init__(self, ia, time_array, mass=None):

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

        """__len__()

        Returns the length of an IonChromatogram.
        """

        return self._ia.size

    def __calc_time_step(self, time_array):

        """__calc_time_step(time_array)

        Calculates the time step.

        @param time_array A numpy object, an array of the retention
        times. 
        """

        td_list = []
        for ii in range(time_array.size-1):
            td = time_array[ii+1]-time_array[ii]
            td_list.append(td) 

        td_array = numpy.array(td_list)
        time_step = td_array.mean()

        return time_step

    def get_intensity_at_index(self, ix):

        """get_intensity_at_index(ix)

        Returns intensity at given index.

        @param ix An integer.
        @return A number.
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > self._ia.size - 1:
            error("index out of bounds")

        return self._ia[ix]

    def get_intensity_array(self):

        """get_intensity_array()

        Returns the entire intensity array.

        @return An array.
        """

        return self._ia

    def get_time_at_index(self, ix):

        """get_time_at_index(ix)

        Returns time at given index.

        @param ix An integer.
        @return A number.
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > self._time_array.size - 1:
            error("index out of bounds")

        return self._time_array[ix]

    def get_time_array(self):

        """get_time_array()

        Returns the entire time array.

        @return An array.
        """

        return self._time_array

    def set_intensity_array(self, ia):

        """set_intensity_array()

        Set the value of intensity array.

        @param ia An array
        @return No value
        """

        self._ia = ia

    def get_time_step(self):

        """get_time_step()

        Returns time step.
        """

        return self._time_step

    def is_tic(self):

        """is_tic()

        Returns True if the ion chromatogram is a total ion
        chromatogram or False otherwise.

        @return A boolean value.
        """

        if self._mass == None:
            return True
        else:
            return False

    def write(self, file_name, minutes=False):

        """write(self, file_name, minutes=False)

        Writes the ion chromatogram to the specified file.

        @param A string representing the file_name to write the ion
            chromatogram to.
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

