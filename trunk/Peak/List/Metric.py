"""
Metric function to compare two peak lists
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

import numpy

from pyms.Peak.Class import Peak

from pyms.Utils.DP import dp
from pyms.Peak.List.DPA.Function import score_matrix_pair
from pyms.Utils.Error import error
from pyms.Utils.Utils import is_list

def metric(peaks1, peaks2, verbose=False, Dw=2.5, Gw=0.30):

    """
    @summary: Compares two peak lists
 
    @param peaks1: A list of Peak objects to compare
    @type peaks1: ListType
    @param peaks2: A list of Peak objects to compare
    @type peaks2: ListType

    @return: metric or distance between two peak lists
    @rtype: FloatType

    @author: Tim Erwin
    """


    if not is_list(peaks1):
        error("First argument of metric() is not a list")

    if not is_list(peaks2):
        error("Second argument of metric() is not a list")

    min_mass,max_mass = get_min_max_mass(peaks1)
    min_mass_exp2,max_mass_exp2 = get_min_max_mass(peaks2)

    if(min_mass_exp2 < min_mass): min_mass = min_mass_exp2
    if(max_mass_exp2 > max_mass): max_mass = max_mass_exp2

    check_mass_spectrum(peaks1,min_mass,max_mass)
    check_mass_spectrum(peaks2,min_mass,max_mass)

    # calculate score matrix for these two runs
    M,B = score_matrix_pair(peaks1, peaks2, Dw)
    # run dp algorithm on this matrix
    result = dp(M, Gw)

    #Calculate Metric
    metric = 0
    match_pos = 0
    gap_penalty=1
    alignment_length = len(result['trace'])
    peak1_pos = 0
    peak2_pos = 0

    metric_display = {}

    for match_value in result['trace']:

      #if match get score from score matrix
      if match_value == 0:
          (peak1,peak2) = result['matches'][match_pos]
	  metric+=M[peak1][peak2]
          rt1 = peaks1[peak1].rt/60.0
          rt2 = peaks2[peak2].rt/60.0
          if(rt1 <= rt2):
              metric_display[rt1] = (rt1,rt2,M[peak1][peak2])
          else:
              metric_display[rt2] = (rt1,rt2,M[peak1][peak2])
	  match_pos += 1
	  peak1_pos += 1
	  peak2_pos += 1
      #else there is no match and apply gap penalty
      else:
          metric+=gap_penalty
          if match_value == 1:
              rt = peaks1[peak1_pos].rt/60.0
              metric_display[rt] = (rt,"-")
              peak1_pos += 1
          if match_value ==2:
              rt = peaks2[peak2_pos].rt/60.0
              peak2_pos += 1
              metric_display[rt] = ("-",rt)

    #Normailise metric
    metric = metric/alignment_length

    if (verbose):
        for key in sorted(metric_display.keys()):
            if (metric_display[key][0] == "-"):
                print "-\t%0.2f" % metric_display[key][1]
            elif(metric_display[key][1] == "-"):
                print "%0.2f\t-" % metric_display[key][0]
            else:
                print "%0.2f\t%0.2f\t%0.2f" % (metric_display[key][0],metric_display[key][1],metric_display[key][2])
        
        print "Metric distance is ",metric

    return metric

def get_min_max_mass(peak_list):

    """
    @summary: Finds the minimum and maximum mass values from a list of peak objects
 
    @param peak_list: A list of Peak objects to compare
    @type peak_list: ListType

    @return: The minimum and Maximum masses
    @rtype: ListType

    @author: Tim Erwin
    """


    min_mass = 1e8
    max_mass = 0
    for peak in peak_list:

        #Check that it is a pyms.Peak.Class.Peak object
	if not (isinstance(peak, Peak)):
            error("Expecting a Peak object")

        #Use mass list
        tmp_min_mass = min(peak.mass_list)
        tmp_max_mass = max(peak.mass_list)
        
        if tmp_min_mass < min_mass: min_mass=tmp_min_mass
        if tmp_max_mass > max_mass: max_mass=tmp_max_mass

    return (min_mass,max_mass)

def check_mass_spectrum(peak_list,min_mass,max_mass):

    """
    @summary: Replaces mass spectrum for each peak in a peak list
              that is bound by minimum and maximum masses
 
    @param peak_list: A list of Peak objects
    @type peak_list: ListType
    @param min_mass: Minimum mass of mass spectrum to create
    @type min_mass: IntType
    @param max_mass: Maximum mass of mass spectrum to create
    @type max_mass: IntType


    @return: none
    @rtype: NoneType

    @author: Tim Erwin
    """

    
    #Create new mass list
    mass_list = []
    for mass in range(min_mass,max_mass+1):
        mass_list.append(mass)

    for peak in peak_list:

        peak_min_mass = peak.mass_list[0]
        
        if(peak.mass_list[0] > min_mass or peak.mass_list[-1]<max_mass):
            peak_offset = peak_min_mass - min_mass

            #Create new mass spectra
            mass_spectrum = numpy.zeros(((max_mass - min_mass)+1), dtype='d')
            mass_spectrum[peak_offset:peak_offset+len(peak.mass_list)] = peak.mass_spectrum[0:]
            peak.mass_spectrum = mass_spectrum
	    peak.mass_list = mass_list

