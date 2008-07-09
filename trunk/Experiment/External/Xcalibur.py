"""
Provides functions for reading exported peaks or peak lists from Xcalibur 
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

import sys, types, os, re, string

from pyms.Peak.Class import Peak
from pyms.Experiment.Class import Experiment
from pyms.IO.ANDI.Class import Xcalibur

from pyms.Utils.IO import file_lines
from pyms.Utils.Error import error

def read_experiments(file_name):

    """ 
    @summary: Reads a file containing Xcalibur experiments into a list of Experiment objects

    @param file_name: The name of the file which lists Xcalibur file names, folled by corresponding netCDF file
    @type file_name: StringType
    @return: A list of Experiment objects
    @rtype: ListType

    @author: Tim Erwin
    @author: Vladimir Likic
    """

    lines = file_lines(file_name)

    experiments = []
    
    for line in lines:

        #Assumes file contains Xcalibur export file followed by corresponding netCDF file
        files = string.split(line," ")
	if(len(files) != 2):
	    error("line must contain Xcalibur file followed by corresponding netCDF file")
        xcalibur_file = string.strip(files[0])
        netcdf_file = string.strip(files[1])
        experiment = load_Xcalibur(xcalibur_file,netcdf_file)

        experiments.append(experiment)

    return experiments

def construct_experiment_from_peaks(file_name):

    """ 
    @summary: Reads a file containing peaks exported from Xcalibur into a Experiment object

    @param file_name: The name of the file which lists contains Xcalibur file names, folled by corresponding netCDF file
    @type file_name: StringType
    @return: A list of Experiment objects
    @rtype: ListType

    @author: Tim Erwin
    @author: Vladimir Likic
    """

    lines = file_lines(file_name)

    peak_list = []

    for peak_file in lines:

        peak_file = string.strip(peak_file)
        peaks = load_peaks(peak_file)
        peak_list = peak_list + peaks

    #Sort peak list
    peak_list.sort(lambda x, y: cmp(x.rt,y.rt))

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment

def load_Xcalibur(file_name, netcdf_file):

    """ 
    @summary: Reads a tab delimited txt file from Xcalibur and creates a pyms.Experiment instance

    @param file_name: The name of the tab delimited txt file exported from Xcalibur
    @type file_name: StringType
    @param netcdf_file: Corresponding netCDF file for exported Xcalibur results
    @type netcdf_file: StringType
    @return: A Experiment object
    @rtype: pyms.Experiment.Class.Experiment

    @author: Tim Erwin
    @author: Vladimir Likic
    """

    lines = file_lines(file_name)
    print " -> Processing Xcalibur file '%s'" % (file_name)
    del lines[0:5] # get rid of the title line

    data = Xcalibur(netcdf_file)
    tic = data.get_tic()

    #Initialise list to store peaks
    peak_list=[]

    for line in lines:

        fields = string.split(line,"\t")

        #initalise variables for peak
        retention_time = float(fields[0])*60.0
        area = fields[3]
        
	#Get mass spectrum at retention time
	index = data.get_index_at_time(retention_time)
        mass_spectrum = data.get_mass_spectrum_at_index(index)
        mass_list = data.get_mass_list()

        #Create peak
        peak = Peak(retention_time, area)
        peak._set_mass_spectrum(mass_spectrum) 
        peak.set_mass_intensity_list(None)
        peak.mass_list = mass_list
        peak_list.append(peak)

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment

def load_peaks(peak_file):

    """ 
    @summary: Loads a file containing one or more peaks exported from Xcalibur

    @param peak_file: The name of the file which contains one or more peaks exported from Xcalibur
    @type peak_file: StringType
    @return: A Peak object
    @rtype: numpy.Peak.Class.Peak

    @author: Tim Erwin
    @author: Vladimir Likic
    """

    #Initialise list to store peaks
    peak_list=[]
    retention_time=0

    lines = file_lines(peak_file)
    print " -> Processing file '%s'" % (peak_file)

    for line in lines:

        #New peak denoted by string beginning with Retention Time
        rt_match = re.search(r'RT: (\d+\.\d*)',line)
        if rt_match: 

            #Store previous peak
            if(retention_time):

                #using dummy area to create peak as we are not using this in
                #the comparison with the metric
                peak = Peak.Class.Peak(retention_time, 1000) 
                peak.set_mass_intensity_list(mass_intensity_list)
                peak_list.append(peak)

            #initalise variables for peak
            retention_time = float(rt_match.group(1))*60.0
            mass_intensity_list={}

        #Get mass list under each peak, intensities are stored as percentages of base peak
        mass_intensity_match = re.search(r'(\d+)\t(\d+)',line)
        if mass_intensity_match:

	    mass = int(mass_intensity_match.group(1))
            intensity = float(mass_intensity_match.group(2))
            mass_intensity_list[mass]=intensity
    
    #Store last peak                
    peak = Peak(retention_time, 1000) 
    peak.set_mass_intensity_list(mass_intensity_list)
    peak_list.append(peak)
    return peak_list


