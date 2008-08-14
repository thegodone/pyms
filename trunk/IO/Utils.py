"""
Utilies for pyms.IO
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
from pyms.IO.Class import IonChromatogram

def is_ionchromatogram(arg):

    """
    @summary: Returns True if the argument is a pyms.IO.Class.IonCromatogram
        object, False otherwise

    @param arg: The argument to be evaluated as IonCromatogram object
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,IonCromatogram):
        return True
    else:
        return False

