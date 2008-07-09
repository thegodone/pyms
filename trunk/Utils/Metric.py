"""
Dynamic Programming routine
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

from pyms.Experiment.Class import Experiment

from pyms.Utils.DP import dp
from pyms.Peak.List.DPA.Function import score_matrix_pair
from pyms.Utils.Error import error

def metric(exp1, exp2, Dw, Gw, sim_method="new"):

    """ Aligns two experiments

    @param A1 The first alignment/experiment
    @param A2 The second alignment/experiment
    @param D "D" parameter
    @param gap Gap parameters
    @param sim_method Use either "new" or "old" methods (default: new)

    @return Merged alignments
    """

    #Error checking arguments
    if not (isinstance(exp1, Experiment)):
        error("First argument is not an experiment object")
    if not (isinstance(exp2, Experiment)):
        error("Second argument is not an experiment object")

    min_mass,max_mass = get_min_max_mass(exp1)
    min_mass_exp2,max_mass_exp2 = get_min_max_mass(exp2)

    if(min_mass_exp2 < min_mass): min_mass = min_mass_exp2
    if(max_mass_exp2 > max_mass): max_mass = max_mass_exp2
    mass_spec(exp1,min_mass,max_mass)
    mass_spec(exp2,min_mass,max_mass)

    # calculate score matrix for these two runs
    M,B = score_matrix_pair(exp1.peaks, exp2.peaks, Dw)
    # run dp algorithm on this matrix
    result = dp(M, Gw)
    #result = dp(M, Gw)

    #Calculate Metric
    metric = 0
    match_pos = 0
    gap_penalty=1
    alignment_length = len(result['trace'])
    exp1_peak_pos = 0
    exp2_peak_pos = 0

    for match_value in result['trace']:
      #if match get score from score matrix
      if match_value == 0:
          (peak1,peak2) = result['matches'][match_pos]
	  metric+=M[peak1][peak2]
          #print match_value, " ",result['matches'][match_pos],M[peak1][peak2]
          rt1 = exp1.peaks[peak1].rt/60.0
          rt2 = exp2.peaks[peak2].rt/60.0
          print "%0.2f\t%0.2f\t%0.2f" % (rt1,rt2,M[peak1][peak2])
	  match_pos += 1
	  exp1_peak_pos += 1
	  exp2_peak_pos += 1
      #else there is no match and apply gap penalty
      else:
          metric+=gap_penalty
          if match_value == 1:
          #if match_value == 2:
              rt = exp1.peaks[exp1_peak_pos].rt/60.0
              print "%0.2f\t-" % rt
              exp1_peak_pos += 1
          if match_value ==2:
          #if match_value ==1:
              rt = exp2.peaks[exp2_peak_pos].rt/60.0
              print "-\t%0.2f" % rt
              exp2_peak_pos += 1

    #Normailise metric
    metric = metric/alignment_length
    print "Metric distance is ",metric

    return

def get_min_max_mass(experiment):

    min_mass = 1e8
    max_mass = 0
    for peak in experiment.peaks:

        if(peak.mass_list):
            #Use mass list
            tmp_min_mass = min(peak.mass_list)
            tmp_max_mass = max(peak.mass_list)
        elif(peak.mass_intensity_list):
            tmp_min_mass = min(peak.mass_intensity_list.keys())
            tmp_max_mass = max(peak.mass_intensity_list.keys()) 
        if tmp_min_mass < min_mass: min_mass=tmp_min_mass
        if tmp_max_mass > max_mass: max_mass=tmp_max_mass

    return min_mass,max_mass

def mass_spec(experiment,min_mass,max_mass):
    
    for peak in experiment.peaks:

        new_mass_spec=[]
        if(peak.mass_intensity_list):
            for i in range(min_mass,max_mass+1):
                if(peak.mass_intensity_list.has_key(i)):
                    new_mass_spec.append(peak.mass_intensity_list[i])
                else:
                    new_mass_spec.append(0.0)
            tmp_mass_spec = numpy.array(new_mass_spec, dtype='d')
            peak._set_mass_spectrum(tmp_mass_spec)
        elif(peak.mass_list):
            if peak.mass_spectrum==None:
                error("Mass spectrum for peak has not been set")
