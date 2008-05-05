"""TopHat.py
 Module TopHat in metab.Baseline
 Provides top-hat baseline corrector. 
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

import sys, copy, math
import numpy
import scipy

from metab.IO import Class
from metab.Utils.Error import error
from metab.Utils.Utils import window_sele_points 
from scipy import ndimage

# default structural element as a fraction of total number of points
_STRUCT_ELM_FRAC = 0.2

def apply(ic, struct=None):

    """apply(ic, struct=None)

    Pefroms top-hat baseline correction.

    @param ic An IonChromatogram object.
    @param struct An integer or string. The width of the structural
         element specified as time string.
    @return Returns a baseline corrected IonChromatogram object.
    """

    if not isinstance(ic, Class.IonChromatogram):
        error("'ic' must be an IonChromatogram object.")
    else:
        ia = copy.deepcopy(ic.get_intensity_array())

    if struct == None:
        struct_pts = int(round(ia.size * _STRUCT_ELM_FRAC))
    else:
        struct_pts = window_sele_points(ic,struct)

    print " -> Top-hat: structural element is %d point(s)" % ( struct_pts )

    str_el = numpy.repeat([1], struct_pts)
    ia = scipy.ndimage.white_tophat(ia, None, str_el)

    ic_bc = copy.deepcopy(ic)
    ic_bc.set_intensity_array(ia)

    return ic_bc 

