"""Utils.py
 Module Utils in pyms.Utils
 Provides general convenience functions.
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

import math, re

import numpy

from Error import error
from types import *

# is_str, is_int, is_float, is_array,  is_number

def is_str(arg):

    """is_str(arg)

    Returns True is argument is string, False otherwise. 
    """

    if type(arg) == str:
        return True 
    else:
        return False

def is_int(arg):

    """is_int(arg)

    Returns True is argument is integer or long, False otherwise. 
    """

    if isinstance(arg,IntType) or isinstance(arg,LongType):
        return True
    else:
        return False

def is_float(arg):

    """is_float(arg)

    Returns True is argument is float, False otherwise. 
    """

    if isinstance(arg,FloatType):
        return True
    else:
        return False

def is_number(arg):

    """is_number(arg)

    Returns True is argument is int, long, or float, False otherwise. 
    """

    if is_int(arg) or is_float(arg):
        return True 
    else:
        return False

def is_list(arg):

    """is_list(arg)

    Returns True is argument is list, typle, or numpy. False
    otherwise. 
    """

    if type(arg) == list or type(arg) == tuple \
            or isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

def is_array(arg):

    """is_array(arg)

    Returns True is argument is numpy, False otherwise. 
    """

    if isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

