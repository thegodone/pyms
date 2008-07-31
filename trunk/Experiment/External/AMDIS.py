"""
Provides functions for reading ELU and FIN files from AMDIS 
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

import sys, types, os, re, string

from pyms.Peak.Class import Peak
from pyms.Experiment.Class import Experiment

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str
from pyms.Utils.IO import file_lines

def read_experiments(file_name):

    """ 
    @summary: Reads a file containing AMDIS experiments into a list of
    Experiment objects

    @param file_name: The name of the file which lists AMDIS file names,
        one file per line
    @type file_name: StringType
    @return: A list of Experiment objects
    @rtype: ListType

    @author: Tim Erwin
    """

    experiment_files = file_lines(file_name)

    experiments = []

    for experiment_file in experiment_files:

        experiment_file = string.strip(experiment_file)
        #TODO read ELU or FIN based on extension, currently expects ELU files
        experiment = load_ELU(experiment_file)
        experiments.append(experiment)

    return experiments

def load_ELU(file_name, uncertain_peaks=0):

    """ 
    @summary: Reads a ELU file from AMDIS and creates a Experiment instance

    @param file_name: A AMDIS ELU file name
    @type file_name: StringType
    @param uncertain_peaks: Boolean flag used to include uncertain masses
        under a peak
    @type uncertain_peaks: BooleanType 
    @return: A  Experiment object
    @rtype: pyms.Experiment.Class.Experiment

    @author: Tim Erwin
    """

    lines = file_lines(file_name)

    print " -> Processing ELU file '%s'" % (file_name)

    #Initialise list to store peaks
    peak_list=[]
    retention_time=0
    area = 0
    
    for line in lines:
    
        #New peak denoted by string beginning with Retention Time
        rt_match = re.search(r'RT(\d+\.\d*)',line)
        if rt_match:

            #Store previous peak
            if(retention_time):

                peak = Peak(retention_time, area) 
                peak.set_mass_intensity_list(mass_intensity_list)
                peak_list.append(peak)

            #initalise variables for peak
            retention_time = rt_match.group(1)
            retention_time=float(retention_time)*60.0
            mass_intensity_list={}

            #Area
            area_match = re.search(r'IS(\d+)',line)
            if area_match: area = area_match.group(1)

            #Abundance of base peak
            bp_match = re.search(r'AM(\d+)',line)
            if bp_match: bp = float(bp_match.group(1))

        #Get mass list under each peak, intensities are stored as percentages
        #of base peak
        mass_intensity_match = re.findall(r'\((\d+)\,(\d+) \)',line)
        for match in mass_intensity_match:
            mass = int(match[0])
            intensity = ( float(match[1])/1000.0 ) * bp
            mass_intensity_list[mass]=intensity


        #Use uncertain peaks
        uncertain_peaks=1
        if uncertain_peaks:
            mass_intensity_ = re.findall(r'(\(\d+\,\d+ B\d\.\d\))',line)
            for match in mass_intensity_match:
                mass = int(match[0])
                intensity = ( float(match[1])/1000.0 ) * bp
                mass_intensity_list[mass]=intensity

    #Create last peak
    peak = Peak(retention_time, area)
    peak.set_mass_intensity_list(mass_intensity_list)
    peak_list.append(peak)

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment

#TODO finish FIN parser can use most of ELU parser code
def load_FIN(file_name, uncertain_peaks=0):

    """ 
    @summary: Reads a FIN file from AMDIS and creates a Experiment instance

    @param file_name: A AMDIS ELU file name
    @type file_name: StringType
    @param uncertain_peaks: Boolean flag used to include uncertain masses
        under a peak
    @type uncertain_peaks: BooleanType 

    @return: An Experiment object
    @rtype: pyms.Experiment.Class.Experiment

    @author: Tim Erwin
    """

    if not is_str(file_name):
        error("'file_name' not a string")
    try:
        fp = open(file_name)
    except IOError:
        error("Cannot open file '%s'" % file_name)

    print " -> Processing FIN file '%s'" % (file_name)


