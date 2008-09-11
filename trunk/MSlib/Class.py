"""
Provides classes related to mass-spectra libraries 
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

class Compound(object):

    """
    @summary: Models a chemical compound

    @author: Saravanan Daylan
    @author: Vladimir Likic
    """

    def __init__(self, cmpd_name):

        """
        @param cmpd_name: The name of the chemical compound
        @type cmpd_name: StringType
        """

        if not is_str(cmpd_name):
            error('Compund name must be a string')

        self.name = cmpd_name
        self.num_peaks = None
        self.mass_spectrum = []
        self.mass_list = []

    def NISTstr2mass(self, mass_spectrum_str):

        """
        @summary: Converts NIST type mass string into mass_list,
        mass_spectrum attributes
        
        @author: Saravanan Daylan
        @author: Vladimir Likic
        """

        if not is_str(mass_spectrum_str):
            error('Argument must be a string')

        tokens = string.split(mass_spectrum_str, ';')
        tokens = tokens[:-1]
        
        if len(tokens) != self.num_peaks:
            error('Error: Number of peaks do not match')
        
        for token in tokens:
            
            pair = string.split(token)
            
            if len(pair) != 2:
                error('Error: Length of pair is != 2')
            
            self.mass_list.append(int(pair[0]))
            self.mass_spectrum.append(int(pair[1]))


