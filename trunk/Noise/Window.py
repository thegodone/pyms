"""Window.py
 Module Window in metab.Noise
 Provides noise smoothing via moving window averaging.
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

import sys
import copy
import numpy

from metab.IO import Class
from metab.Utils.Error import error
from metab.Utils.Utils import window_sele_points 
from metab.Utils.Array import amedian

_DEFAULT_WINDOW = 3

def apply(ic, window=_DEFAULT_WINDOW, median=False):

    """apply(ic, wing_length, median)

    Applies moving window averaging to an IonChromatogram object.

    @param ic An IonChromatogram object.
    @param window An integer or string. Window width selection.
    @param median A boolean. If True median window smoothing will be applied.
    @return An IonChromatogram object
    """

    if not isinstance(ic, Class.IonChromatogram):
        error("'ic' must be an IonChromatogram object.")

    ia = ic.get_intensity_array() # Fetch the intensities.

    wing_length = window_sele_points(ic, window, half_window=True)

    if median:
        ia_denoise = _median_window(ia, wing_length)
    else:
        ia_denoise = _mean_window(ia, wing_length)

    ic_denoise = copy.deepcopy(ic)
    ic_denoise.set_intensity_array(ia_denoise)

    return ic_denoise

def _mean_window(ia, wing_length):

    """_mean_window(ia, wing_length)

    Applies mean-window averaging on the array of intensities.
    @param ia An numpy object. Ion chromatogram intensity array
    @param wing_length An integer value representing the number of
        points on either side of a point in the ion chromatogram
    """

    print " -> Window smoothing (mean): wing is %d point(s)" \
            % (wing_length)

    # The array that will contain the ion chromatogram with noise
    # filtering applied.
    ia_denoise = numpy.repeat([0], ia.size)

    index = 0
    end = ia.size - 1

    while index <= end:
        left = index - wing_length
        right = index + wing_length + 1
        if left < 0: left = 0
        slice = ia[left:right]
        ia_denoise[index] = slice.mean()
        index = index + 1

    return ia_denoise

def _median_window(ia, wing_length):

    """_median_window(ia, wing_length)

    Applies median-window averaging on the array of intensities.
    @param ia An numpy object. Ion chromatogram intensity array
    @param wing_length An integer value representing the number of
        points on either side of a point in the ion chromatogram
    """

    print " -> Window smoothing (median): wing is %d point(s)" \
            % (wing_length)

    # The array that will contain the ion chromatogram with noise
    # filtering applied.
    ia_denoise = numpy.repeat([0], ia.size)

    index = 0
    end = ia.size - 1

    while index <= end:
        left = index - wing_length
        right = index + wing_length + 1
        if left < 0: left = 0
        slice = ia[left:right]
        ia_denoise[index] = amedian(slice)
        index = index + 1

    return ia_denoise

