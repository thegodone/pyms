"""
IO for mass spectral libraries in JCAMP format
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

from Class import MSLibRecord
from pyms.Utils.IO import file_lines
from pyms.Utils.Utils import is_str

def load_jcamp(file_name):

    """
    @summary: Load all records as a list

    @param file_name: The input jcamp file
    @type file_name: StringType

    @return: A records list
    @rtype: ListType

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Reading JCAMP library '%s'" % (file_name)

    lines_list = file_lines(file_name)
    records = []

    #
    # Parse 'lines_list' to produce a list of MSLibRecord objects.
    # This is a mock-up parser that reads the compound name and
    # leaves the mass spectrum blank
    #
    idx1 = 0
    idx2 = 0
    fields = ''
    name_value = ''
    xydata = []
    for line in lines_list:
        prefix = line.find('#')
        if prefix == 0:
            fields = line.split('=')
            if fields[0] == '##CAS NAME':
                idx1 = idx1 + 1
                if idx1 == 1:
                    name_value = fields[1].strip()
                elif idx1 == 2:
                    r = MSLibRecord(name_value, xydata)
                    records.append(r)
                    name_value = fields[1].strip()
                    xydata = []
                    idx1 = 1
                    idx2 = 0
            elif fields[0] == '##XYDATA':
                idx2 = idx2 + 1
        else:
            if prefix == -1 and line != '\n':
                data = line.split()
                data = map(int, data)
                xydata.append(data)
    r = MSLibRecord(name_value, xydata)
    records.append(r)

    return records

