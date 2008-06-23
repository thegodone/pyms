"""
Provides function for reading experiments
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

import string, cPickle

from pyms.Experiment import Class
from pyms.Utils.Utils import *
from pyms.Utils.Error import *

def load_expr(file_name):

    """
    @summary: Loads experiment saved with 'dump'
 
    @param file_name: Experiment file name
    @type file_name: StringType
    @return: none
    @rtype: NoneType
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Loading experiment from the binary file '%s'" % (file_name)

    fp = open(file_name,'r')
    expr = cPickle.load(fp)
    fp.close()

    return expr

def dump_expr(expr, file_name):

    """
    @summary: Dumps an expriment to a file
 
    @param expr: An instance of the Experiment class
    @type expr: Experiment
    @param fle_name: Experiment file name
    @type file_name: StringType
    @return: none
    @rtype: NoneType
    """

    if not isinstance(expr, Class.Experiment):
        error("argument not an instance of the class 'Peak'")

    if not is_str(file_name):
        error("'file_name' not a string")

    fp = open(file_name,'w')
    cPickle.dump(expr, fp, 1)
    fp.close()    

def read_expr(file_name):

    """
    @summary: Reads the experiment files into a list of
        Experiment objects

    @param file_name: The name of the file which lists experiment
        dump file names, one file per line
    @type file_name: StringType
    @return: A list of 'Experiment' instances
    @rtype: ListType
    """

    if not is_str(file_name):
        error("file_name argument must be a string")
    try:
        fp = open(file_name, 'r')
    except IOError:
        error("error opening file '%s' for reading" % file_name)
    
    exprfiles = fp.readlines()
    fp.close()

    exprs = []

    for exprfile in exprfiles:

        exprfile = string.strip(exprfile)
        expr = load_expr(exprfile)

        exprs.append(expr)

    return exprs

