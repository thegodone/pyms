"""IO.py
 Module IO in metab.Utils
 Provides general I/O functions.
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

import types, os, string

from Error import error
from Utils import is_number, is_str, is_list

def open_for_reading(file_name):

    """open_for_reading(file_name)

    Opens file for reading, returns file pointer.

    @param file_name A string.
    """
    if not is_str(file_name):
        error("'file_name' is not a string")
    
    try:
        fp = open(file_name)
    except IOError:
        error("'%s' does not exist" % (file_name))
    return fp

def open_for_writing(file_name):

    """open_for_writing(file_name)

    Opens file for writing, returns file pointer.

    @param file_name A string.
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    try:
        fp = open(file_name, "w")
    except IOError:
        error("Cannot open '%s' for writing" % (file_name))
    return fp

def close_for_reading(fp):

    """close_for_reading(fp)

    Closes file pointer open for reading.

    @param fp A file pointer.
    """

    fp.close()

def close_for_writing(fp):

    """close_for_writing(fp)

    Closes file pointer open for writing.

    @param fp A file pointer.
    """

    fp.close()

def file_lines(file_name):

    """file_lines(file_name)

    Returns lines from a file, as a list.

    @param file_name A string.
    @return A list of lines.
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    fp = open_for_reading(file_name)
    file_lines = fp.readlines()
    close_for_reading(fp)

    return file_lines

def save_data(file_name, data, format_str="%.6f", prepend="", sep=" ",
	compressed=False):

    """save_data(file_name, data, format_str="%.6f", pre="", sep=" ",
	compressed=False)

    Saves a list of numbers or a list of lists of numbers to a file
    with specific formatting.

    @param file_name A string.
    @param data A list of numbers, or a list of lists.
    @param format_str A format string for individual entries.
    @param prepend A string, printed before each row.
    @param sep A string, printed after each number.
    @param compressed A boolean. If True, the output will be gzipped.
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    if not is_list(data):
        error("'data' is not a list")

    if not is_str(prepend):
        error("'prepend' is not a string")

    if not is_str(sep):
        error("'sep' is not a string")

    fp = open_for_writing(file_name)

    # decide whether data is a vector or matrix
    if is_number(data[0]):
        for item in data:
            if not is_number(item):
                error("not all elements of the list are numbers")
        data_is_matrix = 0
    else:
        for item in data:
            if not is_list(item):
                error("not all elements of the list are lists")
        data_is_matrix = 1

    if data_is_matrix:
        for ii in range(len(data)):
            fp.write(prepend)
            for jj in range(len(data[ii])):
                if is_number(data[ii][jj]):
                    fp.write(format_str % (data[ii][jj]))
                    if (jj<(len(data[ii])-1)): fp.write(sep)
                else:
                    error("datum not a number")
            fp.write("\n")
    else:
        for ii in range(len(data)):
            fp.write(prepend)
            fp.write(format_str % (data[ii]))
            fp.write("\n")

    close_for_writing(fp)

    if compressed:
        status = os.system('gzip %s' % (file_name))
        if status != 0:
            error("gzip compress failed")

def read_matrix(file_name, integer=False, compressed=False):

    """read_matrix(file_name, integer=False, compressed=False)

    Reads a matrix from a file and returns a list of lists.
    By default the floats are returned, reading integers
    can be forced with integer=True.

    @param file_name A string.
    @param integer A boolean.
    @param compressed A boolean.
    @return A list of lists, of floats.
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    data_table = read_table(file_name, compressed)

    data_matrix = []
    for row in data_table:
        row_vect = []
        for item in row:
            if integer:
                row_vect.append(int(item))
            else:
                row_vect.append(float(item))
        data_matrix.append(row_vect)

    return data_matrix

def read_table(file_name, compressed=False):

    """read_table(file_name, compressed=False):

    Reads a two dimensional table of strings and returns
    the list of lists.

    @param file_name A string.
    @return A list of lists, of strings.
    """

    if not is_str(file_name):
        error( "file_name not a string")

    data_table = []

    if compressed:
        file_name_gz = file_name + ".gz"
        file_name = ".ReadTable"
        temp_file = file_name + ".gz" 
        status = os.system('cp %s %s' % (file_name_gz, temp_file))
        if (status != 0):
            error( "copy to temporary file failed")
        status = os.system('gzip -d %s' % (temp_file ))
        if (status != 0):
            error("gzip decompress failed")

    lines = file_lines(file_name)

    if compressed:
        status = os.system('rm %s' % (file_name))
        if (status != 0):
            error("removing temporary file failed")

    for line in lines:
        row_vector = []
        items = string.split(line)
        for item in items:
            row_vector.append(item)
        data_table.append(row_vector)

    return data_table

def read_csv_ms(file_name):

    """read_csv_ms(file_name):

    The prototype for a function read CSV matrix saved from
    Microsoft Excel.

    @param file_name A string
    @return A list of lists, containing floats.
    """

    data_table = []

    lines = file_lines(file_name)

    for line in lines:
        row_vector = []
        items = string.split(line,",")
        for item in items:
              s = string.strip(item)
              if s != "":
                  row_vector.append(float(s))
        data_table.append(row_vector)

    return data_table

