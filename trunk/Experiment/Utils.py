"""
Utilities for the subpackage Experiment
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

def cmp_peak_area(peak1, peak2):

    """
    @summary: Compares normalised areas of two peak objects

    @param peak1: A first peak object
    @type peak1: An instance of pyms.Peak.Class.Peak
    @param peak2: A second peak object.
    @type peak2: An instance of pyms.Peak.Class.Peak

    @return: Negative if peak1.norm_area < peak2.norm_area, positive if
         peak1.norm_area > peak2.norm_area, and zero if peak1.norm_area
         == peak2.norm_area
    @rtype: IntType 

    @author: Vladimir Likic
    """

    return cmp(peak1.norm_area, peak2.norm_area)

