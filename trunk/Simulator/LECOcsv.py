"""
Functions to save intensity matrix as text file in LECO csv format
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

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_number, is_str, is_array, is_list, is_int
from pyms.Utils.IO import open_for_writing, close_for_writing

def export_leco_csv(file_name, time_list, mass_list, data):

    """
    @summary: Exports data in LECO CSV format

    @param file_name: File name
    @type file_name: StringType
    @param time_list: List of chromatogram time points
    @type time_list: ListType
    @param mass_list: List of m/z channels
    @type mass_list: ListType
    @param data: Intensity matrix
    @type data: numpy.ndarray

    @return: none
    @rtype: NoneType

    @author: Andrew Isaac
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    if not is_list(time_list):
        error("'time_list' is not a list")

    if not is_list(mass_list):
        error("'mass_list' is not a list")

    if not is_array(data):
        error("'data' is not a numpy ndarray")

    fp = open_for_writing(file_name)

    # Format is text header with:
    # "Scan","Time",...
    # and the rest is "TIC" or m/z as text, i.e. "50","51"...
    # The following lines are:
    # scan_number,time,value,value,...
    # scan_number is an int, rest seem to be fixed format floats.  The format
    # is 0.000000e+000

    # write header
    fp.write("\"Scan\",\"Time\"")
    for ii in mass_list:
        if is_number(ii):
            fp.write(",\"%d\"" % int(ii))
        else:
            error("mass list datum not a number")
    fp.write("\r\n")

    # write lines
    for ii in range(len(time_list)):
        fp.write("%s,%#.6e" % (ii, time_list[ii]))
        for jj in range(len(data[ii])):
            if is_number(data[ii][jj]):
                fp.write(",%#.6e" % (data[ii][jj]))
            else:
                error("datum not a number")
        fp.write("\n")

    close_for_writing(fp)

#    if compressed:
#        status = os.system('gzip %s' % (file_name))
#        if status != 0:
#            error("gzip compress failed")
