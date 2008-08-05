"""
Provides IO functions for peak lists
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

