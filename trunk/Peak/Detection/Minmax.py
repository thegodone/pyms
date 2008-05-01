"""Minmax.py
 Module Minmax in metab.Peak.Detection
 Provides minmax peak detector.
"""

 #############################################################################
 #                                                                           #
 #    META-B software for processing of metabolomic mass-spectrometry data   #
 #    Copyright (C) 2005-6 Vladimir Likic                                    #
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

import sys, copy, math

import numpy
from numpy import linalg

from metab.IO import Class
from metab.Utils.Error import error
from metab.Utils.Utils import *
from metab.Utils.Array import amin

import Integrate
from Constants import *

_DEFAULT_ANGLE_PTS = 3
_DEFAULT_ANGLE_THRESHOLD = 1.0
_DEFAULT_WING_WIDTH = 2
_DEFAULT_PEAK_FACTOR = 10

def apply(ic, noise, window, peak_factor=_DEFAULT_PEAK_FACTOR ):

    """apply(ic, noise, window, peak_factor=_DEFAULT_PEAK_FACTOR)

    Calls minmax peak detection. Returns a list of 'Peak' objects.

    @param ic An IonChromatogram object.
    @param noise A number. The noise estimate. 
    @param window A string. The width of the window for peak
       detection specified as <NUMBER>s or <NUMBER>m, giving a time
       in seconds or minutes, respectively.
    @param peak_factor A number. The value that is multipled to the
        noise level to determine the peak threshold for minimum peak
        intensities.
    @return A list.
    """

    if not isinstance(ic, Class.IonChromatogram):
        error("'ic' must be an IonChromatogram object.")
    else:
        ia = ic.get_intensity_array() # Fetch the intensities.

    if not is_number(noise):
        error("'noise' must be a number.")

    if not is_str(window):
        error("'window' must be a number.")

    if not is_number(peak_factor):
        error("'peak_factor' must be a number.")

    if noise < 0 or peak_factor < 0:
        error("negative argument(s) detected.")

    wing_length = window_sele_points(ic, window, half_window=True)

    print " -> Minmax: using window wing of %d point(s)" %  (wing_length)

    peak_list = _get_peaks(ia, noise, peak_factor, wing_length)
    peaks = Integrate.generate_peaks(ic, peak_list)

    print "    [ found %d peaks ]" % ( len(peaks) )

    return peaks

def _get_peaks(ia, noise, peak_factor, wing_length):

    """_get_peaks(ia, noise, peak_factor, wing_length)

    Peak detection wrapper.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param noise A number. The noise estimate. 
    @param peak_factor A number. The value that is multipled to
        the noise level to determine the peak threshold for
        minimum peak intensities.
    @param wing_length An integer. The width measured in number of
        points to the left and right of a point to be considered
        for an extreme value.
    @return A list of maxima points and their corresponding left and right
    boundaries. Each maximum is modelled using a list containing threex
    elements: the left boundary, maximum point, and right boundary.
    """

    peak_threshold = noise * peak_factor
    maxima_indices = _identify_maxima(ia, wing_length, peak_threshold)
    peak_list = _peak_boundaries(ia, maxima_indices, wing_length)
    _resolve_overlaps(ia, peak_list)
    #_prune_peaks(peak_list, wing_length)
    _adjust_boundaries(ia, peak_list, wing_length, 1.5 ) 

    return peak_list


def _prune_peaks(peak_list, wing_length):

    """_prune_peaks(peak_list, wing_length)

    Rejects peaks for which the difference between the apex and
    left/right boundary is less than wing_length.

    @param peak_list A list of pre-peaks.
    @param wing_length An integer. The width measured in number of
        points to the left and right of a point to be considered
        for an extreme value.
    @return No return value.
    """

    prune_peaks = []

    for peak in peak_list:

        peak_left = peak[LEFT_BOUNDARY]
        peak_apex = peak[MAXIMA]
        peak_right = peak[RIGHT_BOUNDARY]

        if (peak_apex - peak_left) < wing_length or \
           (peak_right - peak_apex) < wing_length:
            prune_peaks.append(peak)

    for peak in prune_peaks:
        peak_list.remove(peak)

def _resolve_overlaps(ia, peak_list):

    """_resolve_overlaps(ia, peak_list)

    Resolves overlapped peaks.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param peak_list A list of pre-peaks.
    @return No return value.
    """

    for ii in range(len(peak_list)-1):

        peak_crnt = peak_list[ii]
        peak_next = peak_list[ii+1]

        if peak_crnt[RIGHT_BOUNDARY] > peak_next[LEFT_BOUNDARY]:
            _adjust_boundary(ia, peak_crnt, peak_next)

def _adjust_boundary(ia, peak_crnt, peak_next):

    """_adjust_boundary(ia, peak_crnt, peak_next)

    Adjusts the boundary between two overlapping peaks. The new
    boundary is set to the minimum between the two peak apexes.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param peak_crnt A list [left,max,right] for peak 1
    @param peak_next A list [left,max,right] for peak 2
    @return No return value.
    """

    crnt_apex = peak_crnt[MAXIMA]
    next_apex = peak_next[MAXIMA]

    minx, minv = amin(ia[crnt_apex+1:next_apex]) 
    split_point = crnt_apex + minx + 1

    peak_crnt[RIGHT_BOUNDARY] = split_point
    peak_next[LEFT_BOUNDARY] = split_point

def _identify_maxima(ia, wing_length, peak_threshold):

    """_identify_maxima(ia, wing_length, peak_threshold)

    Returns an array of indices where each index corresponds to
    a local maximum. A point is a local maximum if the following
    are satisfied:

     1. It is equal or greater than any of 'wing_length' points to
        its left and to its right.
     2. It is greater than at least one of the points to its left,
        and one of the points to its right (out of total 'wing_length'
        points).
     3. A point closer to the edge of 'ia' than 'wing_length' cannot
        be a local maximum.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param wing_length An integer. The width measured in number of
        points to the left and right of a point to be considered
        for an extreme value.
    @param peak_threshold A number. The minimum intensity for a peak,
        maxima below this value are rejected.
    """

    maxima_indices = []

    index = 0
    while index < ia.size:
        is_maxima = _is_local_maximum(ia, index, wing_length)
        is_large_enough = ia[index] >= peak_threshold
        if is_maxima and is_large_enough: # local maximum found 
            maxima_indices.append(index) # append to list
            index = index + wing_length + 1 # index skip for wing_length 
        else:
            index = index + 1

    return maxima_indices


def _peak_boundaries(ia, maxima_indices, wing_length):

    """_peak_boundaries(ia, maxima_indices, wing_length)

    Determines the left and right boundary for each of the maxima.
    Returns a list of maxima points and their corresponding left and
    right boundaries.  Each maxima is modelled using a list of three
    elements: the left boundary, maxium point, and right boundary.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param maxima_indices A numpy object. The list indices of
        the maxima points.
    @param wing_length An integer. The width measured in number of
        points to the left and right of a point to be considered
        for an extreme value.
    @return A list of triplets. Giving indices for [left_bound, maximum,
        right_bound]
    """

    peak_list = []

    # Catch errors in the list of local maxima.
    for crnt_ix in maxima_indices:
        if crnt_ix < wing_length or crnt_ix > ia.size-wing_length-1:
            error("local maxima too close to the boundary")

    for ii in range(len(maxima_indices)-1):
        crnt_ix = maxima_indices[ii]
        next_ix = maxima_indices[ii+1]
        if abs(crnt_ix-next_ix) < wing_length+1:
            error("local maxima too close to one another")

    for ix in range(len(maxima_indices)):

        index = maxima_indices[ix]

        # get the previous and next maximum
        if ix == 0:
            index_prev = 0 
        else:
            index_prev = maxima_indices[ix-1]+1

        if ix == len(maxima_indices)-1:
            index_next = ia.size
        else:
            index_next = maxima_indices[ix+1]

        left_slice = ia[index_prev:index]
        left_slice = left_slice[::-1] # reverse left slice 
        right_slice = ia[index+1:index_next]

        left_min_ix = _get_next_minimum(left_slice, wing_length)
        right_min_ix =  _get_next_minimum(right_slice, wing_length)

        left_boundary = index - left_min_ix - 1
        right_boundary = index + right_min_ix + 1

        peak_list.append([left_boundary, index, right_boundary])

    return peak_list

def _get_next_minimum(ia_slice, wing_length):

    """_get_next_minimum(ia_slice, wing_length)

    Determines the first next minimum in 'ia_slice' based on 'wing_length'.
    'ia_slice' is processed as index increases until the first minimum
    is found. If no minimum is found, the minimum is set to the first
    point after 'wing_length' counting from the end of the slice.

    @param ia_slice An numpy object. Ion chromatogram intensity
        array slice.
    @param wing_length An integer. The width measured in number of
        points to the left and right of a point to be considered
        for an extreme value.
    @return An integer. The index of the minimum value.
    """

    boundary_index = None

    for index in range(ia_slice.size):

        if _is_local_minimum(ia_slice, index, wing_length):
            boundary_index = index
            break

    # if the minimum was not found, set the boundary to the last
    # point in the slice. If the two peaks are exactly 'wing_length'
    # points appart (the closest possible), the right boundary of
    # first peak will be set to the left adjacent point of the apex
    # of the second peak. And the left boundary of the second peak
    # will be set to the right adjacent point of the apex of the
    # first peak apex. If a peak is at the beginning/end of a signal,
    # and no local minima are found when searching towards the boundary
    # the boundary will be set to the first/last point of the signal.
    if boundary_index == None:
        boundary_index = ia_slice.size - 1

    return boundary_index

##
# Returns True if the point is a local maximum or False otherwise. 
# A point is a local maximum if the following are satisfied:
#    1. It is equal or greater than any of 'wing_length' points to
#       its left and to its right.
#    2. It is greater than at least one of the points to its left,
#       and one of the points to its right (out of total 'wing_length'
#       points).
#    3. A point closer to the edge of 'ia' than 'wing_length' cannot
#       be a local maximum.
#
# @param ia An numpy object, ion chromatogram intensity array.
# @param wing_length An integer. The width measured in number of
#     points to the left and right of a point to be considered
#     for extreme value.
# @return A boolean value.

def _is_local_maximum(ia, index, wing_length):

    """_is_local_maximum(ia, index, wing_length)

    Returns True if the point is a local maximum or False otherwise. 
    A point is a local maximum if the following are satisfied:
       1. It is equal or greater than any of 'wing_length' points to
          its left and to its right.
       2. It is greater than at least one of the points to its left,
          and one of the points to its right (out of total 'wing_length'
          points).
       3. A point closer to the edge of 'ia' than 'wing_length' cannot
          be a local maximum.

    @param ia An numpy object. Ion chromatogram intensity array 
        slice.
    @param index An integer. Index at which to evaluate if local maximum. 
    @return A boolean value.
    """

    if index > ia.size:
        error("The index (%d) is out of bounds (max index: %d)." % \
              (index, ia.size - 1))

    intensity = ia[index]

    # Get the index for left/right boundary, fold the edges.
    left_index = index - wing_length
    if left_index < 0: left_index = 0
    left_slice = ia[left_index:index]
    
    right_index = index + wing_length + 1
    if right_index > ia.size: right_index = ia.size
    right_slice = ia[index+1:right_index]

    # Test if local maximum. Points closer to the edges than
    # 'wing_length' cannot be a maximum.
    if left_slice.size < wing_length:
        left_wing_valid = False
    else:
        left_wing_valid = _is_greater_than(left_slice, intensity)

    if right_slice.size < wing_length:
        right_wing_valid = False 
    else:
        right_wing_valid = _is_greater_than(right_slice, intensity)

    # Return value: True if both wings passed 
    if left_wing_valid and right_wing_valid:
        return True
    else:
        return False

def _is_greater_than(ia, intensity):

    """_is_greater_than(ia, intensity)

    Returns True if 'intensity' is greater than or equal to every
    element in 'ia' _and_ is greater than at least one element.

    @param ia An numpy object. Ion chromatogram intensity array 
        slice.
    @param intensity A number. Intensity to which slice is compared. 
    """

    if ia.size < 1:
        error("empty array for comparison")

    # numpy.where(ia <= intensity)[0] returns an array of indices
    # refering to elements of 'ia' that are smaller or equal than
    # 'intensity'
    lteq = numpy.where(ia <= intensity)[0]
    lt = numpy.where(ia < intensity)[0]

    if lteq.size == ia.size and lt.size > 0:
        return True
    else:
        return False

def _is_local_minimum(ia, index, wing_length):

    """_is_local_minimum(ia, index, wing_length)

    Returns True if the point is a local minimum or False otherwise. 
    A point is a local maximum if the following are satisfied:
       1. It is equal or smaller than any of 'wing_length' points to
          its left and to its right.
       2. It is smllaer than at least one of the points to its left,
          and one of the points to its right (out of total 'wing_length'
          points).
       3. A point closer to the edge of 'ia' than 'wing_length' cannot
          be a local minimum.

    @param ia An numpy object. Ion chromatogram intensity array 
        slice.
    @param index An integer. Index at which to evaluate if local maximum. 
    @return A boolean value.
    """

    if index > ia.size:
        error("The index (%d) is out of bounds (max index: %d)." % \
              (index, ia.size - 1))

    intensity = ia[index]

    # Get the index for left/right boundary, fold the edges.
    left_index = index - wing_length
    if left_index < 0: left_index = 0
    left_slice = ia[left_index:index]
    
    right_index = index + wing_length + 1
    if right_index > ia.size: right_index = ia.size
    right_slice = ia[index+1:right_index]

    # Test if local maximum. Points closer to the edges than
    # 'wing_length' cannot be a maximum.
    if left_slice.size < wing_length:
        left_wing_valid = False
    else:
        left_wing_valid = _is_less_than(left_slice, intensity)

    if right_slice.size < wing_length:
        right_wing_valid = False
    else:
        right_wing_valid = _is_less_than(right_slice, intensity)

    # Return value: True if both wings passed 
    if left_wing_valid and right_wing_valid:
        return True
    else:
        return False

def _is_less_than(ia, intensity):

    """_is_less_than(ia, intensity)

    Returns True if 'intensity' is less than or equal to every
    element in 'ia' _and_ is less than at least one element.

    @param ia An numpy object. Ion chromatogram intensity array 
        slice.
    @param intensity A number. Intensity to which slice is compared. 
    """

    if ia.size < 1:
        error("empty array for comparison")

    # numpy.where(ia >= intensity)[0] returns an array of indices
    # refering to elements of 'ia' that are greater or equal to 
    # 'intensity'
    lteq = numpy.where(ia >= intensity)[0]
    lt = numpy.where(ia > intensity)[0]

    if lteq.size == ia.size and lt.size > 0:
        return True
    else:
        return False

def _adjust_boundaries(ia, peak_list, angle_width=_DEFAULT_ANGLE_PTS, \
        angle_threshold=_DEFAULT_ANGLE_THRESHOLD):

    """ _adjust_boundaries(ia, peak_list, angle_width, angle_threshold)

    Adjusts peak boundaries by fitting the straight line through the
    edge points, and trimming point for which the line formes an angle
    with the time axis which is less than the threshold.

    @param ia An numpy object. Ion chromatogram intensity array.
    @param peak_list A list of pre-peaks.
    @param angle_width An integer. The number of point for the best fit. 
    @param angle_threshold A number. The angle threshold for peak boundary
        points.  A line of best fit drawn through an arbitrary number of
        points at the peak boundaries must form an angle greater than
        the angle threshold.  If not, the boundary point will be trimmed.
    @return No return value. 
    """

    for peak_item in peak_list:

        left_index = peak_item[LEFT_BOUNDARY]
        apex_index = peak_item[MAXIMA]
        right_index = peak_item[RIGHT_BOUNDARY]

        apex = ia[apex_index]
        left_slice = ia[left_index:apex_index]
        right_slice = ia[apex_index+1:right_index+1]
        # Reverse the numpy. _trim_boundary() expects the
        # numpy from boundary to apex.
        right_slice = right_slice[::-1]

        left_slice = _trim_boundary(left_slice, apex, angle_width, \
                angle_threshold)
        right_slice = _trim_boundary(right_slice, apex, angle_width, \
                angle_threshold)

        peak_item[LEFT_BOUNDARY] = apex_index - left_slice.size
        peak_item[RIGHT_BOUNDARY] = apex_index + right_slice.size

def _trim_boundary(ia_slice, apex, angle_width, angle_threshold):

    """_trim_boundary(ia_slice, apex, angle_width, angle_threshold)

    Performs peak boundary trimming. 

    @param ia_slice An numpy object. Ion chromatogram intensity
        array slice representing initial peak boundary.
    @param apex A number. Intensity at peak apex .
     @param angle_width An integer. The number of point for the best fit. 
    @param angle_threshold A number. The angle threshold for peak boundary
        points.  A line of best fit drawn through an arbitrary number of
        points at the peak boundaries must form an angle greater than
        the angle threshold.  If not, the boundary point will be trimmed.
    @return A numpy object. Trimmed peak boundary.
    """

    if ia_slice.size <= angle_width: return ia_slice

    index = 0

    while index + angle_width < ia_slice.size:

        x1 = index
        x2 = index + angle_width
        x_values = numpy.arange(x1, x2+1)
        y_values = ia_slice[x_values]

        gradient = _best_fit_line(y_values, x_values)[0]
        rise = angle_width * gradient
        rise = rise / apex
        run = angle_width
        angle = math.atan(rise/run)
        angle = math.degrees(angle)
        index = index + 1

        if angle > angle_threshold: break

    return ia_slice[index-1:]


def _best_fit_line(y_values, x_values=None):

    """_best_fit_line(y_values, x_values=None)

    Calculates the coefficients to a linear line of best fit.  The

    @param y_values A numpy object representing the y-values 
    @param x_values A numpy object representing the x-values
    @return A tuple. Best fit coefficients where the first element is
        the gradient, and the second is the y-intercept.
    """

    if x_values == None: x_values = numpy.arange(y_values.size)

    x_array = numpy.repeat([1], x_values.size)
    x_array = numpy.array([x_values, x_array])
    x_array = numpy.transpose(x_array)
    y_array = copy.deepcopy(y_values)

    alpha = numpy.dot(numpy.transpose(x_array), x_array)
    beta = numpy.dot(numpy.transpose(x_array), y_array)
    ialpha = linalg.inv(alpha)
    coefficients = numpy.dot(ialpha, beta)

    return coefficients[0], coefficients[1]

