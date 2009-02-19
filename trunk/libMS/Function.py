"""
General functions for library matching 
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

import numpy

from pyms.Utils.Error import error

def ms_lib_match(ms_lib, ms):

    """
    @summary: Matches mass spectrum 'ms' against all the records
        in the mass spectral library 'ms_lib'

    @para ms: The mass spectrum object
    @type ms: pyms.IO.Class.MassSpectrum

    @return: The compound name, similarity value, and mass spectrum
    @rtype: ListType

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    # for each record, set the attribute 'mass_spectrum'
    max = 0.0
    compound = ''
    mass = []
    mass_list = ms.mass_list
    im = ms.mass_spec
    for record in ms_lib.records:
        score = 0.0
        if len(record.mass_record) > 0:
            array1 = im
            array2 = []
            for i in range(len(array1)):
                array2.append(0)
            for item in record.mass_record:
                if item[0] >= mass_list[0] and item[0] <= mass_list[len(mass_list) - 1]:
                    array2[item[0] - 50] = item[1]
            mass_spect1 = numpy.array(array1, dtype='d')
            mass_spect2 = numpy.array(array2, dtype='d')
            mass_spect1_sum = numpy.sum(mass_spect1 ** 2, axis=0)
            mass_spect2_sum = numpy.sum(mass_spect2 ** 2, axis=0)
            top = numpy.dot(mass_spect1, mass_spect2)
            bot = numpy.sqrt(mass_spect1_sum * mass_spect2_sum)
            cos = top / bot
            score = cos
        if score > max:
            max = score
            compound = record.name
            mass = record.mass_record
    result = [compound, max, mass]
    return result
