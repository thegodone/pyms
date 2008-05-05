"""Analysis.py
 Module Window in pyms.Noise
 Noise estimator based on randomly placed windows.
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

import sys, random, copy, math

import numpy

from pyms.IO import Class 
from pyms.Utils.Error import error
from pyms.Utils.Array import amad
from pyms.Utils.Utils import window_sele_points 

_DEFAULT_WINDOW = 256
_DEFAULT_N_WINDOWS = 1024

def apply(ic, window=_DEFAULT_WINDOW, n_windows=_DEFAULT_N_WINDOWS, rand_seed=None ):

    """apply(ic, window, n_windows, rand_seed)

    Applies noise analysis based on randomly placed windows.

    The noise value is calculated by repeatedly and picking
    random windows (of a specified width) and calculating
    median absolute deviation (MAD).  The noise estimate is
    given by the minimum MAD.

    @param ic An IonChromatogram object.
    @param window An integer or string. Window width selection.
    @param n_windows An integer. The number of windows to calculate
    @param rand_seed A number. The seed value used in the random
        number generator will be multiplied by this number.
    @return A floating point value. The noise estimate
    """

    if not isinstance(ic, Class.IonChromatogram):
        error("'ic' must be an IonChromatogram object.")

    ia = ic.get_intensity_array() # Fetch the intensities.
    
    # Create an instance of the Random class
    if rand_seed != None:
        generator = random.Random(rand_seed)
    else:
        generator = random.Random()

    window_pts = window_sele_points(ic,window)

    maxi = ia.size - window_pts
    est_noise = math.fabs(ia.max()-ia.min()) 
    best_window_pos = None
    seen_positions = []

    cntr = 0
    
    while cntr < n_windows:
        # generator.randrange(): last point not included in range
        try_pos = generator.randrange(0, maxi+1)
        # Only process the randomly generated window if the window
        # is unique.
        if try_pos not in seen_positions:
            end_slice = try_pos + window_pts
            slice = ia[try_pos:end_slice]
            crnt_mad = amad(slice)
            if crnt_mad < est_noise:
                est_noise = crnt_mad
                best_window_pos = try_pos
        cntr = cntr + 1
        seen_positions.append(try_pos)
  
    return est_noise

