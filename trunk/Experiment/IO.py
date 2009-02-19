"""
Functions related to experiment input/output
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

from pyms.IO.ANDI.Class import Xcalibur, ChemStation
from pyms.Peak.List.IO import read_xcalibur_peaks, read_amdis_peaks, read_analyzerpro_peaks
from pyms.Utils.Error import error 
from pyms.Experiment.Class import Experiment 
from pyms.Utils.Utils import is_str
from pyms.Utils.IO import file_lines

def load_expr(file_name):

    """
    @summary: Loads an experiment saved with 'dump_expr'
 
    @param file_name: Experiment file name
    @type file_name: StringType

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
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
    @summary: Dumps expriment to a file
 
    @param expr: An experiment object
    @type expr: pyms.Experiment.Class.Experiment
    @param file_name: The name of the file
    @type file_name: StringType

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
    """

    if not isinstance(expr, Experiment):
        error("argument not an instance of the class 'Experiment'")

    if not is_str(file_name):
        error("'file_name' not a string")

    fp = open(file_name,'w')
    cPickle.dump(expr, fp, 1)
    fp.close()

    print " -> Experiment '%s' saved as '%s'" % (expr.expr_code, file_name)

def read_expr_list(file_name):

    """
    @summary: Reads the set of experiment files and returns a list of
    Experiment objects

    @param file_name: The name of the file which lists experiment
        dump file names, one file per line
    @type file_name: StringType

    @return: A list of Experiment instances
    @rtype: ListType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("file_name argument must be a string")
    try:
        fp = open(file_name, 'r')
    except IOError:
        error("error opening file '%s' for reading" % file_name)
    
    exprfiles = fp.readlines()
    fp.close()

    exprl = []

    for exprfile in exprfiles:

        exprfile = string.strip(exprfile)
        expr = load_expr(exprfile)

        exprl.append(expr)

    return exprl

def load_xcalibur_expr(file_name, netcdf_file):

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
    print " -> Processing Xcalibur experiment"
    
    peak_list = read_xcalibur_peaks(file_name)
    andi_data = Xcalibur(netcdf_file)

    for peak in peak_list:
        peak.set_mass_spectrum(andi_data)

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment

def read_xcalibur_expr_list(file_name):

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
        xcalibur_peak_file = string.strip(files[0])
        netcdf_file = string.strip(files[1])
        experiment = load_xcalibur_expr(xcalibur_peak_file,netcdf_file)
        experiments.append(experiment)

    return experiments


def load_amdis_expr(file_name, andi_data=None):

    """ 
    @summary: Reads a ELU file from AMDIS and creates a pyms.Experiment instance

    @param file_name: The name of a ELU file from AMDIS
    @type file_name: StringType
    @param andi_data: Optional ANDI-MS data to use for mass spectrum instead of AMDIS mass data
    @type andi_data: IO.ANDI.Class.ANDIMS_reader object
    @return: A Experiment object
    @rtype: pyms.Experiment.Class.Experiment

    @author: Tim Erwin
    @author: Vladimir Likic
    """
    print " -> Processing AMDIS experiment"
    
    peak_list = read_amdis_peaks(file_name)

    if(andi_data):
        for peak in peak_list:
            peak.set_mass_spectrum(andi_data)

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment

def load_analyzerpro_expr(file_name, netcdf_file):

    """ 
    @summary: Loads a peak list exported from AnalyzerPro

    @param file_name: The name of the txt file exported from AnalyzerPro
    @type file_name: StringType
    @param netcdf_file: Corresponding netCDF file for the AnalyzerPro results
    @type netcdf_file: StringType
    @return: A Experiment object
    @rtype: pyms.Experiment.Class.Experiment

    @author: Tim Erwin
    @author: Vladimir Likic
    """
    print " -> Processing AnalyzerPro experiment"

    peak_list = read_analyzerpro_peaks(file_name)
    andi_data = ChemStation(netcdf_file)

    for peak in peak_list:
        peak.set_mass_spectrum(andi_data)

    #Create Experiment object
    experiment = Experiment(file_name, peak_list)

    return experiment
