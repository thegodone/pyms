"""
Classes used for the manipulation of mass spectral libraries in JCAMP format
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

import IO

class MSLib(object):

    """
    @summary: Models a MS library

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    def __init__(self, file_name):

        """
        @summary: Initialize the library

        @para file_name: The input jcamp file name
        @type file_name: StringType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        self.records = IO.load_jcamp(file_name)
        self.mass_list = None

    def set_for_search(self, mass_list):

        """
        @summary: Sets the mass spectrum value for each record

        @para mass_list: The mass list values along the m/z
        @type mass_list: ListType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        self.mass_list = mass_list

        # for each record, set the attribute 'mass_spectrum'
        for record in self.records:
            record._set_mass_spectrum(mass_list)

class MSLibRecord(object):

    """
    @summary: Models a MS libarary Record

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    def __init__(self, name, mass_record):

        """
        @summary: Initialize the record

        @para name: The compound name
        @type name: StringType
        @para mass_record: The list mass spectrum values
        @type mass_record: ListType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        self.name = name
        self.mass_record = mass_record
        self.mass_spectrum = []

    def _set_mass_spectrum(self, mass_list):

        """
        @summary: Sets mass spectrum value for each compound

        @para mass_list: The mass list values along the m/z
        @type mass_list: ListType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        #
        # set the value of record.mass_spectrum here
        #
        for item in self.mass_record:
            if item[0] >= mass_list[0] and item[0] <= mass_list[len(mass_list) - 1]:
                self.mass_spectrum.append(item[1])
            else:
                self.mass_spectrum.append(0)