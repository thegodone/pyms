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

    # initialise best match attributes
    best_score = 0.0
    compound = None
    mass = []

    # find min/max mz in the query mass spectrum
    min_mz = min(ms.mass_list)
    max_mz = max(ms.mass_list)

    # fill in the the mz_index dictionary. Keys are mz values,
    # values are indices in the intensity list
    mz_index = {}
    for ii in range(len(ms.mass_list)):
        mz = ms.mass_list[ii]
        mz_index[mz] = ii

    # loop over all records in the library
    for record in ms_lib.records:

        record_score = 0.0

        # bomb-out if record has an empty mass spectrum
        if len(record.mass_record) == 0:
            error("library record with empty mass spectrum")

        # initialize record mass spectrum
        record_mass_spec = []
        for i in range(len(ms.mass_list)):
            record_mass_spec.append(0)

        # fill in record mass spectrum
        for item in record.mass_record:

            record_mz = item[0]
            intensity = item[1]

            # consider only record mz values that are in the range
            # of the query mass spectrum
            if (record_mz >= min_mz) and (record_mz <= max_mz):

                # query mass spectrum must have this mz
                if not record_mz in mz_index.keys():
                    error("mz value missing in the library record")
                else: # assign intensity to the record mz value
                    ii = mz_index[record_mz]
                    record_mass_spec[ii] = intensity

        # get the mass spectral similarity score
        record_score = costheta(ms.mass_spec, record_mass_spec)

        # consider if the best hit
        if record_score > best_score:
            best_score = record_score
            compound = record.name
            mass = record.mass_record

    result = [compound, best_score, mass]

    return result

def costheta(v1, v2):

    """
    """

    mass_spect1 = numpy.array(v1, dtype='d')
    mass_spect2 = numpy.array(v2, dtype='d')
    mass_spect1_sum = numpy.sum(mass_spect1 ** 2, axis=0)
    mass_spect2_sum = numpy.sum(mass_spect2 ** 2, axis=0)

    top = numpy.dot(mass_spect1, mass_spect2)
    bot = numpy.sqrt(mass_spect1_sum * mass_spect2_sum)
    cos = top / bot

    return cos


