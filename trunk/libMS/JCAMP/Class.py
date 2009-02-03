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
import numpy
from time import time

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

    def printl(self, begin=0, end=None):

        """
        """

        if end == None:
            end = len(self.records)

        for ii in range(begin,end):
            record = self.records[ii]
            print "(%d)" % ( ii ),
            print record.name
            for item in record.mass_record:
                print "\t", item

    def match(self, mass_list, im):

        """
        @summary: Matches a scan point against all the records

        @para mass_list: The mass list values along the m/z
        @type mass_list: ListType
        @para im: The mass spectrum of the scan point
        @type im: ListType

        @return: The compound name and the similarity value
        @rtype: StringType, FloatType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        # for each record, set the attribute 'mass_spectrum'
        max = 0.0
        compound = ''
        for record in self.records:
            score = 0.0
            if len(record.mass_record) > 0:
                array1 = []
                array2 = []
                for item in record.mass_record:
                    if item[0] >= mass_list[0] and item[0] <= mass_list[len(mass_list) - 1]:
                        array1.append(item[1])
                        array2.append(im[mass_list.index(item[0])])
                mass_spect1 = numpy.array(array1, dtype='d')
                mass_spect2 = numpy.array(array2, dtype='d')
                mass_spect1_sum = numpy.sum(mass_spect1 ** 2, axis=0)
                mass_spect2_sum = numpy.sum(mass_spect2 ** 2, axis=0)
                top = numpy.dot(mass_spect1, mass_spect2)
                bot = numpy.sqrt(mass_spect1_sum * mass_spect2_sum)
                cos = top / bot
                score = 1.0 - cos
            if score > max:
                max = score
                compound = record.name
        return compound, max

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
