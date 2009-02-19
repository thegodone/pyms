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
from pyms.Utils.Utils import is_str, is_int
from pyms.Peak.Class import Peak
from pyms.Utils.IO import open_for_writing, close_for_writing, file_lines

from Utils import is_peak_list

_VERBOSE_LEVEL = "verbose_level"
_DEFINED_VERBOSE_LEVELS = [0, 1]
_TIME_UNITS = "time_units"
_MINUTES = "minutes"
_SECONDS = "seconds"

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
    @summary: Reads AMDIS ELU file, and returns the list of peak objects

    @param file_name: AMDIS .ELU file name
    @type file_name: StringType
    @param uncertain_masses: Use uncertain masses in mass spectra
    @type file_name: BoolType

    @return: Returns the list of pyms.Peak.Class.Peak objects
    @rtype: ListType

    @author: Tim Erwin
    """

    lines = file_lines(file_name)

    print " -> Reading AMDIS ELU file '%s'" % (file_name)

    #Initialise list to store peaks
    peak_list=[]
    peaks={}
    unique_peaks={}
    retention_time=0
    area = 0
    new_peak_name = None
    new_peak_percentage = None    

    for line in lines:
    
        #New peak denoted by string beginning with Retention Time
        rt_match = re.search(r'RT(\d+\.\d*)',line)
        if rt_match:

            peak_name = new_peak_name
            peak_percentage = new_peak_percentage

            #Get peak name and percentage match
            #name_match = re.search(r'SC(\d+)',line)
            percentage_match = re.search(r'%(\d+\.\d*)',line)

            #if not (name_match and percentage_match):
            if not (percentage_match):
                error("ELU is not of expected format")
            else:
                #new_peak_name = name_match.group(1)
                new_peak_percentage = percentage_match.group(1)

            #Store previous peak only if best match
            #there may be multiple models of a single peak
            if (retention_time):

                peak = Peak(retention_time, area)
                mass_list,mass_spectrum = create_mass_spec(mass_intensity_list)
                peak.mass_list = mass_list
                peak.mass_spectrum = mass_spectrum

                if peaks.has_key(peak_name):
                    #Store best peak
                    if(peak_percentage > peaks[peak_name]):
                        unique_peaks[peak_name] = peak
                        peaks[peak_name] = peak_percentage
                else:
                    unique_peaks[peak_name] = peak
                    peaks[peak_name] = peak_percentage

            #initalise variables for peak
            retention_time = rt_match.group(1)
            retention_time=float(retention_time)*60.0
            new_peak_name = retention_time
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
    peak = Peak(retention_time, area)
    mass_list,mass_spectrum = create_mass_spec(mass_intensity_list)
    peak.mass_list = mass_list
    peak.mass_spectrum = mass_spectrum

    if peaks.has_key(peak_name):
        #Store best peak
        if(peak_percentage > peaks[peak_name]):
            unique_peaks[peak_name] = peak
            peaks[peak_name] = peak_percentage
    else:
        #print "Storing Peak ",retention_time
        unique_peaks[peak_name] = peak
        peaks[peak_name] = peak_percentage

    return unique_peaks.values()

def read_analyzerpro_peaks(file_name):

    """
    @summary: Reads Xcalibur peak report, and returns the list of peak objects

    @param file_name: Peak list file name
    @type file_name: StringType

    @return: Returns the list of pyms.Peak.Class.Peak objects
    @rtype: ListType

    @author: Tim Erwin
    """

    lines = file_lines(file_name)
    print " -> Reading AnalyzerPro peak file '%s'" % (file_name)
    del lines[0:3] # get rid of the title line

    #Initialise list to store peaks
    peak_list=[]

    for line in lines:

        fields = string.split(line," ")

        #initalise variables for peak
        retention_time = float(fields[1])*60.0
        peak_area = float(fields[5])

        #Create peak
        peak = Peak(retention_time, peak_area)
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

def write_peaks(peaks, file_name, minutes=False, verbose_level=0):

    """
    Writes a peak list to a file.
 
    @param peaks A list. List of Peak objects.
    @param fle_name A string. Peak list file name.
    @param minutes A boolean. If True, convert time to minutes. 
    @param verbose_level An integer. Verbose level.
    @return Returns no values. 
    """

    if not is_peak_list(peaks):
        error("'peaks' not a list of peak objects")

    if not is_str(file_name):
        error("'file_name' not a string")

    if not is_int(verbose_level):
        error("'verbose_level' must be an integer.")
    else:
        if not verbose_level in _DEFINED_VERBOSE_LEVELS:
            error("unknown 'verbose_level'.")

    level_string = "# %s %d" % (_VERBOSE_LEVEL, verbose_level)

    if minutes:
        time_string = "# %s %s" % (_TIME_UNITS, _MINUTES)
    else:
        time_string = "# %s %s" % (_TIME_UNITS, _SECONDS)

    print " -> Writing peak list to '%s' (verbose level %d)" % \
        (file_name, verbose_level)

    fp = open_for_writing(file_name)

    fp.write("# PyMS peak list file\n")
    fp.write("%s\n" % (level_string))
    fp.write("%s\n" % (time_string))

    if verbose_level == 0:

        fp.write("# peak_num, rt_apex, raw_area\n");

        cntr = 0

        for peak in peaks:
            cntr = cntr + 1
            rt_apex = peak.rt
            raw_area = peak.raw_area
            if minutes:
               rt_apex = rt_apex/60.0
            fp.write("%4d %10.3f %12.1f\n" % (cntr, rt_apex, raw_area))

    elif verbose_level == 1:

        fp.write("# peak_num, pt_left, pt_apex, pt_right,");
        fp.write(" rt_apex, raw_area\n");

        cntr = 0

        for peak in peaks:

            cntr = cntr + 1

            if peak.pt_bounds == None:
               print "\t peak pt_bounds not defined."
               error("cannot write peak list at verbose level %d" % \
                   (verbose_level))
            else:
                pt_left = peak.pt_bounds[0]
                pt_apex = peak.pt_bounds[1]
                pt_right = peak.pt_bounds[2]

            rt_apex = peak.rt

            if minutes:
               rt_apex = rt_apex/60.0

            raw_area = peak.raw_area

            fp.write("%4d %5d %5d %5d %10.3f %12.1f\n" % (cntr, \
                    pt_left, pt_apex, pt_right, rt_apex, raw_area))

    close_for_writing(fp)

