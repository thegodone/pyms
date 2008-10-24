"""
Provides input/output functions for peak lists
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

import string
import re
import numpy

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str 
from pyms.Utils.IO import file_lines
from pyms.Peak.Class import Peak

def read_chem_station_peaks(file_name):

    """
    @summary: Reads ChemStation peak integration report, and returns
    the list of peak objects

    @param file_name: Peak list file name
    @type file_name: StringType

    @return: Returns the list of pyms.Peak.Class.Peak objects
    @rtype: ListType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Reading ChemStation peak integration report '%s'" % (file_name)

    lines_list = file_lines(file_name)
    lines_list = lines_list[:len(lines_list)-5] # Cut the comment

    lines_to_remove = []
    for line in lines_list:
        fields = string.split(line)
        if len(fields) == 0:
            lines_to_remove.append(line)
    for line in lines_to_remove:
        lines_list.remove(line)

    peaks = []
    read = 0
    peak_id = -1
    cntr = 0

    for line in lines_list:

        cntr = cntr + 1
        fields = string.split(line)

        if (len(fields)>1) and (fields[0] == "---"):
            read = 1

        if len(fields)==0:
            read = 0

        if read:
            peak_id = peak_id + 1

            if peak_id:

                fields2 = string.split(line)
                peak_rt = float(fields2[1])

                # -2 correction was found empirically to
                # match the peak apex
                pt_left = int(fields2[2])-2
                pt_apex = int(fields2[3])-2
                pt_right = int(fields2[4])-2
                pt_bounds = [pt_left, pt_apex, pt_right]

                short_line = line[33:]

                fields2 = string.split(short_line)
                raw_area = int(fields2[1])

                if len(fields2) == 5:
                    peak_tag = string.lower(fields2[4])
                else:
                    peak_tag = None

                peak = Peak(peak_rt, raw_area, minutes=True)

                # set the peak tag and peak boundaries in points
                peak.set_peak_tag(peak_tag)
                peak.set_pt_bounds(pt_bounds)

                peaks.append(peak)

    return peaks

def read_xcalibur_peaks(file_name):

    """
    @summary: Reads Xcalibur peak report, and returns the list of peak objects

    @param file_name: Peak list file name
    @type file_name: StringType

    @return: Returns the list of pyms.Peak.Class.Peak objects
    @rtype: ListType

    @author: Tim Erwin
    """

    lines = file_lines(file_name)

    print " -> Reading Xcalibur peak file '%s'" % (file_name)
    del lines[0:5] # get rid of the title line

    #Initialise list to store peaks
    peak_list=[]

    for line in lines:

        fields = string.split(line,"\t")

        #initalise variables for peak
        retention_time = float(fields[0])*60.0
        peak_area = float(fields[3])

        #Create peak
        peak = Peak(retention_time, peak_area)
        peak_list.append(peak)

    return peak_list

def read_amdis_peaks(file_name,uncertain_masses=False):

    """
    @summary: Reads Xcalibur peak report, and returns the list of peak objects

    @param file_name: Peak list file name
    @type file_name: StringType

    @return: Returns the list of pyms.Peak.Class.Peak objects
    @rtype: ListType

    @author: Tim Erwin
    """

    lines = file_lines(file_name)

    print " -> Reading AMDIS ELU file '%s'" % (file_name)

    #Initialise list to store peaks
    peak_list=[]
    retention_time=0
    area = 0
    peak_name = None
    peak_percentage = None    

    for line in lines:
    
        #New peak denoted by string beginning with Retention Time
        rt_match = re.search(r'RT(\d+\.\d*)',line)
        if rt_match:

            prev_peak_name = peak_name
            prev_peak_percentage = peak_percentage

            #Get peak name and percentage match
            name_match = re.search(r'SC(\d+)',line)
            percentage_match = re.search(r'%(\d+\.\d*)',line)

            if not (name_match and percentage_match):
                error("ELU is not of expected format")
            else:
                peak_name = name_match.group(1)
                peak_percentage = percentage_match.group(1)

            #Store previous peak only if best match (there may be multiple models of
            #a single peak)
            if (retention_time):
                if (prev_peak_name == peak_name) and \
                        (prev_peak_percentage < peak_percentage):
                    peak = Peak(retention_time, area) 
                    mass_list,mass_spectrum = create_mass_spec(mass_intensity_list)
                    peak.mass_list = mass_list
                    peak.mass_spectrum = mass_spectrum
                    peak_list.append(peak)

            #initalise variables for peak
            retention_time = rt_match.group(1)
            retention_time=float(retention_time)*60.0
            mass_intensity_list={}

            #Area
            area_match = re.search(r'IS(\d+)',line)
            if area_match: area = float(area_match.group(1))

            #Abundance of base peak
            bp_match = re.search(r'AM(\d+)',line)
            if bp_match: bp = float(bp_match.group(1))

        #Get mass list under each peak, intensities are stored as percentages of base peak
        mass_intensity_match = re.findall(r'\((\d+)\,(\d+) \)',line)
        for match in mass_intensity_match:
            mass = int(match[0])
            intensity = ( float(match[1])/1000.0 ) * bp
            mass_intensity_list[mass]=intensity

        #Use uncertain peaks
        if uncertain_masses:
            mass_intensity_ = re.findall(r'(\(\d+\,\d+ B\d\.\d\))',line)
            for match in mass_intensity_match:
                mass = int(match[0])
                intensity = ( float(match[1])/1000.0 ) * bp
                mass_intensity_list[mass]=intensity

    #Create last peak
    if(prev_peak_name == peak_name) and (prev_peak_percentage < peak_percentage):
        peak = Peak(retention_time, area)
        mass_list,mass_spectrum = create_mass_spec(mass_intensity_list)
        peak.mass_list = mass_list
        peak.mass_spectrum = mass_spectrum
        peak_list.append(peak)

    return peak_list

def create_mass_spec(mass_intensity_list):

    mass_values = mass_intensity_list.keys()
    min_mass = min(mass_values)
    max_mass = max(mass_values)

    mass_spectrum = numpy.zeros(((max_mass - min_mass)+1), dtype='d')
    mass_list = []
    i = 0

    for mass in range(min_mass,max_mass+1):
        mass_list.append(mass)
        if (mass_intensity_list.has_key(mass)):
            mass_spectrum[i] = mass_intensity_list[mass]
        i+=1

    return mass_list, mass_spectrum

