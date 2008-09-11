"""
Provides classes related to NIST libraries 
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
from pyms.Utils.IO import file_lines
from pyms.MSlib.Class import Compound

class NIST_library(object):

    """
    @summary: NIST library reader 

    @author Saravanan Daylan
    @author Vladimir Likic
    """ 

    def __init__(self, file_name):

        """
        @param file_name: The name of the NIST library file
        @type file_name: StringType
        """

        if not is_str(file_name):
            error('File name must be a string')
 
        self.compounds = self.parse_NIST(file_name)

    def parse_NIST(self, file_name):

        """
        @summary: Parses NIST library file

        @param file_name: The name of the NIST library file
        @type file_name: StringType

        @author Saravanan Daylan
        @author Vladimir Likic
        """ 

        __CMPD_NAME_KEYWORD = "NAME"
        __NUM_PEAKS_KEYWORD = "NUM PEAKS"

        nist_lines = file_lines(file_name)

        compounds = []
        crnt_cmpd = None 
        prev_cmpd = None
        new_cmpd_flag = False
        collect_ms_flag = False
        mass_spectrum_string = ""
        first_mass_line_flag = False
        cmpd_counter = 0

        for line in nist_lines:

            fields = string.split(line, ":")

            if len(fields)>0 and fields[0].upper() == __CMPD_NAME_KEYWORD:

                prev_cmpd = crnt_cmpd

                keyword_value = fields[1]
                cmpd_name = keyword_value.strip()
                crnt_cmpd = Compound(cmpd_name)
                compounds.append(crnt_cmpd)
                cmpd_counter =+ 1
                
                new_cmpd_flag = True
                collect_ms_flag = False

            elif len(fields)>0 and fields[0].upper() == __NUM_PEAKS_KEYWORD:

                keyword_value = fields[1]
                num_peaks_str = keyword_value.strip()
                crnt_cmpd.num_peaks = int(num_peaks_str)

                collect_ms_flag = True
                first_mass_line_flag = True
                
            if prev_cmpd != None and new_cmpd_flag:
                prev_cmpd.NISTstr2mass(mass_spectrum_string)
                mass_spectrum_string = ""
                new_cmpd_flag = False

            if collect_ms_flag:
                if not first_mass_line_flag:
                    mass_spectrum_string = mass_spectrum_string + line.strip()
                else:
                    first_mass_line_flag = False
            
        if cmpd_counter > 0:
            crnt_cmpd.NISTstr2mass(mass_spectrum_string)
        
        return compounds 

