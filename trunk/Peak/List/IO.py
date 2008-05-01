"""IO.py
 Module IO in metab.Peak.List
 Provides IO functions for peak lists.
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

import os, string, cPickle

from metab.Peak import Class
from metab.Utils.Error import error
from metab.Utils.IO import file_lines, open_for_writing, close_for_writing
from metab.Peak.List.Utils import is_peak_list 
from metab.Utils.Utils import *
from metab.IO.ANDI.Class import ChemStation

_VERBOSE_LEVEL = "verbose_level"
_DEFINED_VERBOSE_LEVELS = [0, 1, 7]

_TIME_UNITS = "time_units"
_MINUTES = "minutes"
_SECONDS = "seconds"

def read(file_name):

    """read(file_name)

    Reads meta-b peak list.

    @param fle_name A string. Peak list file name.
    @return Returns list of Peak objects.
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Reading peak list '%s'" % (file_name),

    lines_list = file_lines(file_name)

    # Determine verbose_level.
    fields = string.split(lines_list[1])

    if fields[1] != _VERBOSE_LEVEL:
        error("cannot find verbose level tag")
    else:
        verbose_level = int(fields[2])

        if not verbose_level in _DEFINED_VERBOSE_LEVELS:
            error("unknown verbose level (%d)" % (verbose_level))

    print "(verbose level %s)" % (verbose_level)

    # Determine time units.
    fields = string.split(lines_list[2])

    if fields[1] != _TIME_UNITS:
        error("cannot find time untis tag")
    else:
        if fields[2] == _MINUTES:
            minutes = True
        elif fields[2] == _SECONDS:
            minutes = False
        else:
            error("unknown time units specification ('%s')" % (fields[2]))

    # Cut the header.
    lines_list = lines_list[4:]

    peaks = []

    for line in lines_list:

        fields = string.split(line)

        if verbose_level == 0:

            rt_apex = float(fields[1])
            raw_area = float(fields[2])

            if minutes: rt_apex = rt_apex*60.0

            peak = Class.Peak(rt_apex, raw_area)

            peaks.append(peak)

        elif verbose_level == 1:

            pt_left = int(fields[1])
            pt_apex = int(fields[2])
            pt_right = int(fields[3])
            pt_bounds = [pt_left, pt_apex, pt_right]

            rt_apex = float(fields[4])

            raw_area = float(fields[5])

            if minutes: rt_apex = rt_apex*60.0

            peak = Class.Peak(rt_apex, raw_area)
            peak.set_pt_bounds(pt_bounds)

            peaks.append(peak)

        elif verbose_level == 7:

            pt_left = int(fields[1])
            pt_apex = int(fields[2])
            pt_right = int(fields[3])
            pt_bounds = [pt_left, pt_apex, pt_right]

            rt_left = float(fields[4])
            rt_apex = float(fields[5])
            rt_right = float(fields[6])
           
            intensity = float(fields[7])
            raw_area = float(fields[8])

            if minutes:
                rt_left = rt_left*60.0
                rt_apex = rt_apex*60.0
                rt_right = rt_right*60.0

            rt_bounds = [rt_left, rt_apex, rt_right]

            peak = Class.Peak(rt_apex, raw_area)

            peak.set_intensity(intensity)
            peak.set_pt_bounds(pt_bounds)
            peak.set_rt_bounds(rt_bounds)

            peaks.append(peak)

        else:
            error("this statement should never execute")

    return peaks

def load(file_name):

    """load(file_name)

    Loads peaks saved with 'dump'.
 
    @param fle_name A string. Peak list file name.
    @return Returns no values.
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Loading peaks from the binary file '%s'" % (file_name)

    fp = open(file_name,'r')
    peaks = cPickle.load(fp)
    fp.close()

    return peaks

def dump(peaks, file_name):

    """dump(peaks, file_name)

    Dumps a peak list to a file.
 
    @param peaks A list. List of Peak objects.
    @param fle_name A string. Peak list file name.
    @return Returns no values. 
    """

    if not is_peak_list(peaks):
        error("'peaks' not a list of peak objects")

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Saving peaks binary file '%s'" % (file_name)
    
    fp = open(file_name,'w')
    cPickle.dump(peaks, fp, 1)
    fp.close()

def write(peaks, file_name, minutes=False, verbose_level=0):

    """write(peaks, file_name, minutes=False, verbose_level=0)

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

    fp.write("# meta-b peak list file\n")
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

    elif verbose_level == 7:

        fp.write("# peak_num, pt_left, pt_apex, pt_right, rt_left,");
        fp.write(" rt_apex, rt_right, intensity, raw_area\n");

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

            if peak.rt_bounds == None:
                print "\t peak rt_bounds not defined."
                error("cannot write peak list at verbose level %d" % \
                        (verbose_level))
            else:
                rt_left = peak.rt_bounds[0]
                rt_apex = peak.rt
                rt_right = peak.rt_bounds[2]

            raw_area = peak.raw_area

            if peak.intensity == None:
                print "\t peak intensity not defined."
                error("cannot write peak list at verbose level %d" % \
                        (verbose_level))
            else: 
                intensity = peak.intensity

            if minutes:
               rt_left = rt_left/60.0
               rt_apex = rt_apex/60.0
               rt_right = rt_right/60.0

            fp.write("%4d %5d %5d %5d %10.3f %10.3f %10.3f %10d %12.1f\n" % \
                   (cntr, pt_left, pt_apex, pt_right, rt_left, rt_apex,
                   rt_right, intensity, raw_area))

    close_for_writing(fp)

def read_chem_station(file_name):

    """read_chem_station(file_name)

    Reads ChemStation peak integration report, and returns the list
    of peak objects.

    @param file_name A string. Peak list file name.
    @return Returns the list of Peak objects.
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Reading ChemStation peak integration report '%s'" % (file_name)

    lines_list = file_lines(file_name)
    lines_list = lines_list[:len(lines_list)-5] # Cut the comments.

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
                peak_rt = float(fields2[1])*60.0

                # -2 correction was found empirically to
                # match the peak apex
                pt_left = int(fields2[2])-2
                pt_apex = int(fields2[3])-2
                pt_right = int(fields2[4])-2
                pt_bounds = [pt_left, pt_apex, pt_right]

                short_line = line[33:]

                fields2 = string.split(short_line)
                peak_raw_area = int(fields2[1])

                if len(fields2) == 5:
                    peak_tag = string.lower(fields2[4])
                else:
                    peak_tag = None

                peak = Class.Peak(peak_rt, peak_raw_area)
                peak.set_peak_tag(peak_tag)
                peak.set_pt_bounds(pt_bounds)

                peaks.append(peak)

    return peaks

