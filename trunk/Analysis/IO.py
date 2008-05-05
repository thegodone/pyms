"""IO.py
 Module IO in pyms.Analysis
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

from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import * 
from pyms.Utils import IO
from pyms.Utils import Array 

def split_two_groups(data_list, groups):

    """split_two_groups(data_list, groups):

    Splits the data in 'data_list' into two groups, according to
    labels given in the list 'groups'.

    @param data_list A list.
    @param grups A list.
    @return (data_list1, data_list2) A tuple.
    """

    data_dict = {}

    for ii in range(len(data_list)):
        group = groups[ii]
        if not data_dict.has_key(group):
           data_dict[group] = []
        data_dict[group].append(data_list[ii])

    groups = data_dict.keys()

    if len(groups) != 2:
        error("the number of groups not equal 2")

    data_list1 = []

    for item in data_dict[groups[0]]:
        data_list1.append(item)

    data_list2 = []

    for item in data_dict[groups[1]]:
        data_list2.append(item)

    return data_list1, data_list2

def load_data_groups(annot_file):

    """load_data_groups(annot_file)

    Reads data groups from the annotation file.

    @param annot_file A string
    @return group A list
    """

    lines = IO.file_lines(annot_file)

    cntr = 0
    label_str = None
    label_str_prev = None
    label_crnt = 0
    groups = []

    for line in lines:
        fields = string.split(line,"/")
        label_str = fields[1]
        if label_str != label_str_prev:
            label_crnt = label_crnt + 1
            label_str_prev = label_str
        groups.append(label_crnt)

    return groups

def load_ptype_cvs(cvs_file_name):

    """load_ptype_cvs(cvs_file_name)

    Extracts the matrix of peak areas and the matrix of retention
    times from the DPA prototype output.

    @param cvs_file_name A string
    @return No return value
    """

    if not is_str(cvs_file_name):
        error("file name not a string")
    
    lines = IO.file_lines(cvs_file_name)
    tmp = lines.pop(0) # get rid of the title line

    data_matrix = []
    rt_matrix = []

    for line in lines:

        data_row = []
        rt_row = []

        fields = string.split(line,",")

        if len(fields) % 2 == 0: 
            error("even number of columns in CSV file")

        tmp = fields.pop(0) # get rid of row number

        for ii in range(len(fields)):
            if ii % 2 == 0:
                if fields[ii].strip() == "--":
                    value = 0.0
                else:
                    value = float(fields[ii])
                rt_row.append(value)
            else:
                if fields[ii].strip() == "--":
                    value = 0.0
                else:
                    value = float(fields[ii])
                data_row.append(value)

        data_matrix.append(data_row)
        rt_matrix.append(rt_row)

    return data_matrix, rt_matrix 


