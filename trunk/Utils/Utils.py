"""Utils.py
 Module Utils in metab.Utils
 Provides general convenience functions.
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

import math, re

import numpy

from Error import error
from types import *

def is_str(arg):

    """is_str(arg)

    Returns True is argument is string, False otherwise. 
    """

    if type(arg) == str:
        return True 
    else:
        return False

def is_int(arg):

    """is_int(arg)

    Returns True is argument is integer or long, False otherwise. 
    """

    if isinstance(arg,IntType) or isinstance(arg,LongType):
        return True
    else:
        return False


def is_float(arg):

    """is_float(arg)

    Returns True is argument is float, False otherwise. 
    """

    if isinstance(arg,FloatType):
        return True
    else:
        return False

def is_number(arg):

    """is_number(arg)

    Returns True is argument is int, long, or float, False otherwise. 
    """

    if is_int(arg) or is_float(arg):
        return True 
    else:
        return False

def is_list(arg):

    """is_list(arg)

    Returns True is argument is list, typle, or numpy. False
    otherwise. 
    """

    if type(arg) == list or type(arg) == tuple \
            or isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

def is_array(arg):

    """is_array(arg)

    Returns True is argument is numpy, False otherwise. 
    """

    if isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

def list_or_array(arg):

    """list_or_array(arg)

    Returns True is argument is list or numpy, False otherwise. 
    """

    if is_list(arg) or is_array(arg):
        return True 
    else:
        return False

def is_string_number(arg):

    """is_string_number(arg)

    This function determines if its argument, x, is in the format of a
    number. It can be number can be in integer, floating point,
    scientific, or engineering format. The function returns True if
    the argument is formattted like a number, and False otherwise.
    (by gyro funch from Active State Python Cookbook.)
    """

    NUM_RE = re.compile(r'^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$')

    if NUM_RE.match(str(arg)):
        return True
    else:
        return False

def window_sele_points(ic, window_sele, half_window=False):

    """window_sele_points(ic, window_sele, half_window)

    Converts 'window_sele' to points based on the times step in IC.

    @param ic An IonChromatogram object.
    @param window_sele An integer or a string. The window specification. 
        if integer, taken as the number of points. If a string, must 
        of the form "<NUMBER>s" or "<NUMBER>m", specifying a time
        in seconds or minutes, respectively.
    @param half_window  A boolean. If True, half-window (window wing)
       is returned.
    """

    if not is_int(window_sele) and not is_str(window_sele):
        error("'window' must be an integer or a string")

    if is_int(window_sele):

        if half_window:
            if window_sele % 2 == 0: 
                error("window must be an odd number of points")
            else: 
                points = int(math.floor(window_sele*0.5))
        else:
            points = window_sele
    else:
        time = time_str_seconds(window_sele)

        time_step = ic.get_time_step()

        if half_window: time = time*0.5

        points = int(math.floor(time/time_step))

    if half_window:
        if points < 1: error("window too small (half window=%d)" % (points))
    else:
        if points < 2: error("window too small (window=%d)" % (points))

    return points

def time_str_seconds(time_str):

    """time_str_seconds(time_str)

    Resolves time string of the form "<NUMBER>s" or "<NUMBER>m",
    returns time in seconds.

    @param time_str A string. Must be of the form "<NUMBER>s" or
        "<NUMBER>m" where "<NUMBER>" is a valid number.
    @return A number. Time in seconds.
    """

    if not is_str(time_str):
        error("time string not a string")

    time_number = time_str[:-1]
    time_spec = time_str[-1].lower()

    if not is_string_number(time_number):
       print " --> received time string '%s'" % (time_number)
       error("improper time string")

    if not time_spec == "s" and not time_spec == "m":
        error("time string must end with either 's' or 'm'")

    time = float(time_number)

    if time_spec == "m":
        time = time*60.0
    
    return time

def time_secs_to_mins(secs):

    """time_secs_to_mins(secs)

    Converts numeric seconds to minutes string, to .1s precision. For example,
    381.32s = '6:21.3'

    @param secs Seconds, numeric argument.
    @return A string. Time in seconds.
    """

    if not is_number(secs):
        error("seconds is not a numeric argument")

    return '%d:%02d.%d' % (secs / 60, secs % 60, secs % 60 % 1.0 * 10)

