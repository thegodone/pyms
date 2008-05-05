"""Utils.py
 Module Utils in pyms.Peak.List
 Utilities for peak lists.
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

import math

from pyms.Peak import Class
from pyms.Utils.Error import error
from pyms.Utils.Utils import *

def is_peak_list(peaks):

    """is_peak_list(peaks)

    Returns True is 'peaks' is valid peak list.

    @param peaks A list of peak objects
    @return A boolean
    """

    flag = True

    if not is_list(peaks):
        flag = False
    else:
        for item in peaks:
           if not isinstance(item, Class.Peak):
              flag = False

    return flag

def find_peak_by_rt(rt, peaks, tol_str):

    """find_peak_by_rt(rt, peaks, tol_str)

    Finds the closest peak to given retention time. The closest
    peak must:
        1. Be closer to the target retention time than tolerance.
        2. Be the only peak which is closer to the target retention
           time than tolerance.

    @param rt A number. Target retention time in seconds.
    @param peaks A list. The list of candidate peaks.
    @param tol_str A time string, "<NUMBER>s" or "<NUMBER>m"
        to specify a cutoff in seconds or minutes.
    @return A peak object
    """

    if not is_number(rt):
        error("'rt' not a number")

    if not is_peak_list(peaks):
        error("'peaks' not a peak list")

    if not is_str(tol_str):
        error("'tol_str' not a string")

    matches = []

    cutoff = time_str_seconds(tol_str)

    for peak in peaks:

        rt_diff = math.fabs(rt - peak.rt)

        if rt_diff < cutoff:
            matches.append(peak)

    if len(matches) == 0 or len(matches) > 1:
        peak_match = None
    else:
        peak_match = matches[0]

    return peak_match

def sele_peaks_by_rt(peaks, rt_range):

    """sele_peaks_by_rt(peaks, rt_range)

    Selects peaks by retention time range.

    @param peaks A list of peak objects
    @param rt_range A list of two time strings, specifying lower and
           upper retention time
    @return A list of peak objects
    """

    if not is_peak_list(peaks):
        error("'peaks' not a peak list")

    if not is_list(rt_range):
        error("'rt_range' not a list")
    else:
        if len(rt_range) != 2:
            error("'rt_range' must have exactly two elements")

        if not is_str(rt_range[0]) or not is_str(rt_range[1]):
            error("lower/upper retention time limits must be strings")

    rt_lo = time_str_seconds(rt_range[0])
    rt_hi = time_str_seconds(rt_range[1])

    if not rt_lo < rt_hi:
        error("lower retention time limit must be less than upper")

    peaks_sele = []

    for peak in peaks:
        if peak.rt > rt_lo and peak.rt < rt_hi:
            peaks_sele.append(peak)

    print " -> %d peaks selected" % (len(peaks_sele))

    return peaks_sele

def cmp_by_rt(peak1, peak2):

    """cmp_by_rt(peak1, peak2)

    Compares two peak objects by retention time.

    @param peak1 A peak object
    @param peak2 A peak object
    @return An integer Negative, zero, or positive depending on how
            peak1 and peak2 retention times compare (built-in cmp
            function)
    """

    return cmp(peak1.rt, peak2.rt)

