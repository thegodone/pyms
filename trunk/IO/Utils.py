"""
Utilies for pyms.IO
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

from pyms.Utils.Error import error
from pyms.IO.Class import IonChromatogram
from pyms.Utils.Utils import is_int, is_str
from pyms.Utils.Time import time_str_secs

def is_ionchromatogram(arg):

    """
    @summary: Returns True if the argument is a pyms.IO.Class.IonCromatogram
        object, False otherwise

    @param arg: The argument to be evaluated as IonCromatogram object
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,IonChromatogram):
        return True
    else:
        return False

def ic_window_points(ic, window_sele, half_window=False):

    """
    @summary: Converts window selection parameter into points based on
        the time step in an ion chromatogram 

    @param ic: ion chromatogram object relevant for the conversion
    @type ic: pyms.IO.Class.IonChromatogram
    @param window_sele: The window selection. If integer, taken as the
        number of points. If a string, must of the form "<NUMBER>s" or
        "<NUMBER>m", specifying a time in seconds or minutes, respectively
    @type window_sele: IntType or StringType 
    @param half_window: Specifies whether to return half-window
    @type half_window: Booleantype

    @author: Vladimir Likic
    """

    if not is_int(window_sele) and not is_str(window_sele):
        error("'window' must be either an integer or a string")

    if is_int(window_sele):

        if half_window:
            if window_sele % 2 == 0: 
                error("window must be an odd number of points")
            else: 
                points = int(math.floor(window_sele*0.5))
        else:
            points = window_sele
    else:
        time = time_str_secs(window_sele)
        time_step = ic.get_time_step()

        if half_window:
            time = time*0.5

        points = int(math.floor(time/time_step))

    if half_window:
        if points < 1: error("window too small (half window=%d)" % (points))
    else:
        if points < 2: error("window too small (window=%d)" % (points))

    return points

