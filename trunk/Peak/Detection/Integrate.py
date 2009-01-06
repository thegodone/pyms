"""
Functions common to all peak integrators
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

from pyms.Peak.Class import Peak 

from Constants import *

def generate_peaks(ic, peak_list):

    """
    @summary: Generates a list of peak objects from the list of pre-peaks,
        and integrates peak areas

    @param ic: An IonChromatogram object
    @type ic: numpy.ndarray
    @param peak_list: A list of pre-peaks
    @type peak_list: ListType

    @return: A list of peak objects
    @rtype: ListType

    -----------------------------------------------------------

    generate_peaks(ic, peak_list)

    Generates a list of peak objects from the list of pre-peaks,
    and integrates peak areas.

    @param ic An IonChromatogram object.
    @param peak_list A list of pre-peaks.
    @return A list of peak objects
    """

    ia = ic.get_intensity_array()

    peaks = []

    for peak_item in peak_list:

        pt_left = peak_item[LEFT_BOUNDARY]
        pt_apex = peak_item[PEAK_APEX]
        pt_right = peak_item[RIGHT_BOUNDARY]

        # peak areas is the sum of intensity from peak left
        # to its right boundary (boundary points included) 
        area_slice = ia[pt_left:pt_right+1]
        raw_area = area_slice.sum()

        rt_apex = ic.get_time_at_index(pt_apex)
        peak = Peak(rt_apex, raw_area)
        peak.set_pt_bounds([pt_left, pt_apex, pt_right])

        peaks.append(peak)

    return peaks