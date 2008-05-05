"""Integrate.py
 Module Integrate in metab.Peak.Detection
 Provides functions common to peak detectors.
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

from metab.Peak import Class
from Constants import *

def generate_peaks(ic, peak_list):

    """generate_peaks(ic, peak_list)

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
        pt_apex = peak_item[MAXIMA]
        pt_right = peak_item[RIGHT_BOUNDARY]
        pt_bounds = [pt_left, pt_apex, pt_right]

        intensity = ic.get_intensity_at_index(pt_apex)

        rt_left = ic.get_time_at_index(pt_left)
        rt_apex = ic.get_time_at_index(pt_apex)
        rt_right = ic.get_time_at_index(pt_right)
        rt_bounds = [rt_left, rt_apex, rt_right]

        # peak areas is the sum of intensity from peak left
        # to its right boundary (boundary points included) 
        area_slice = ia[pt_left:pt_right+1]
        raw_area = area_slice.sum()

        peak = Class.Peak(rt_apex, raw_area)
        peak.set_intensity(intensity)
        peak.set_pt_bounds(pt_bounds)
        peak.set_rt_bounds(rt_bounds)

        peaks.append(peak)

    return peaks

