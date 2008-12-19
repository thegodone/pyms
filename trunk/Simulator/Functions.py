"""
GCMS_simulator utility functions
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

import copy
import random
import numpy
import math

"""
Provides a function simulate pure mass spectrum given the number of m/z
channels num_channels
"""
def mass_spectrum(num_channels):

    """
    Returns a simulated pure mass spectrum given the number
    of m/z channels num_channels

    A mass spectrum contains between 1 and N/10 components
    (the precise number is randomly chosen), with random
    intensity in range 0-1.
    """

    #% initialize the m/z vector
    M = numpy.zeros((num_channels), 'd')

    #% the number of non-zero m/z values
    # NB: P < num_channels
    P = int((num_channels/10)*random.random());

    #% generate mass spectrum
    rp = range(num_channels)
    random.shuffle(rp)
    for ii in range(P):
        kk = rp[ii]
        M[kk] = random.random()
#        M[kk] = 1

    return M

"""
Provides a function to return the gaussian for a given point
"""
def gaussian(point, mean, width, scale):
    return scale*math.exp(-((point-mean)**2)/(2*width**2))

"""
Provides a function to simulate ion chromatogram of a pure component
"""
def chromatogram(num_chrom_points):

    """
    Returns a simulated ion chromatogram of a pure component
    The ion chromatogram contains a single gaussian peak.
    The peak position is randomly chosen, but not to be closer
    to the edges for more than the fraction 'offset' of the
    total time domain.
    The peak width is randomly chosen between 'w1' and 'w2'.
    The peak intensity is randomly chosen between 0 and 1.
    """

    #% peak width limits
    w1 = 0.0005
    w2 = 0.0015

    #% peak position offset from the edges
    offset = 0.02
    d = offset
    g = 1-offset

    V = numpy.zeros((num_chrom_points), 'd')

    peak_width = w1 + (w2-w1)*random.random()
    peak_pos = d + (g-d)*random.random()
    peak_scale = random.random()

    for ii in range(num_chrom_points):
        x = float(ii)/num_chrom_points    # float to avoid int division
        V[ii] = gaussian(x, peak_pos, peak_width, peak_scale)

    return V
