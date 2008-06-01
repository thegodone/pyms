"""Utils.py
 Module Utils in pyms.Utils
 Provides general utility functions.
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

import types

import numpy

from Error import error

def is_str(arg):

    """
    Returns True if the argument is a string, False otherwise. 

    @type arg: arbitrary 
    @param arg: The argument to be evaluated as a string. 
    @trype: BooleanType
    @return: A boolean indicator True or False.
    """

    if isinstance(arg,types.StringType):
        return True 
    else:
        return False

def is_int(arg):

    """
    Returns True if the argument is an integer, False otherwise. 

    @type arg: arbitrary 
    @param arg: The argument to be evaluated as an integer.
    @trype: BooleanType 
    @return: A boolean indicator True or False.
    """

    if isinstance(arg,types.IntType) or isinstance(arg,types.LongType):
        return True
    else:
        return False

def is_float(arg):

    """
    Returns True if the argument is a float, False otherwise. 

    @type arg: arbitrary 
    @param arg: The argument to be evaluated as a float.
    @trype: BooleanType 
    @return: A boolean indicator True or False.
    """

    if isinstance(arg,types.FloatType):
        return True
    else:
        return False

def is_number(arg):

    """
    Returns True if the argument is a number (integer or float),
    False otherwise. 
   
    @type arg: arbitrary
    @param arg: The argument to be evaluated as a number.
    @trype: BooleanType 
    @return: A boolean indicator True or False.
    """

    if is_int(arg) or is_float(arg):
        return True 
    else:
        return False

def is_list(arg):

    """
    Returns True if the argument is a list, tuple, or numpy array,
    False otherwise.

    @type arg: arbitrary
    @param arg: The argument to be evaluated as a list.
    @trype: BooleanType 
    @return: A boolean indicator True or False.
    """

    if isinstance(arg,types.ListType) or isinstance(arg,types.TupleType) \
            or isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

def is_array(arg):

    """
    Returns True if the argument is a numpy array, False otherwise. 

    @type arg: arbitrary 
    @param arg: The argument to be evaluated as a numpy array. 
    @trype: BooleanType 
    @return: A boolean indicator True or False.
    """

    if isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

