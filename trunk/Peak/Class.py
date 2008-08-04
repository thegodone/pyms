"""
Contains Peak-related classes
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

import math

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_int, is_number, is_list

class Peak:

    """
    @summary: Models a signal peak

    A peak object is initialised with retention time and raw peak area.

    @author: Vladimir Likic
    """

    def __init__(self, peak_rt=0.0, raw_area=0, minutes=False):

        """
        @param peak_rt: Retention time
        @type peak_rt: A number
        @param raw_area: Raw peak area
        @type raw_area: A number
        @param minutes: Retention time units flag. If True, retention time
            is in minutes; if Flase retention time is in seconds
        @type minutes: BooleanType 
        """

        if not is_number(peak_rt):
            error("'peak_rt' must be a number")

        if not is_number(raw_area):
            error("'raw_area' must be a number")

        if minutes:
            peak_rt = peak_rt/60.0

        self.rt = float(peak_rt)
        self.raw_area = float(raw_area)

        # optional attributes
        self.mass_spectrum = None 
        self.mass_list = None 
        self.tag = None

    def set_pt_bounds(self, pt_bounds):

        """
        @summary: Sets peak boundaries in points

        @param pt_bounds: A list containing left, apex, and right
            peak boundary (points)
        @type pt_bounds: ListType

        @return: none
        @rtype: NoneType
        """

        if not is_list(pt_bounds):
            error("'pt_bounds' must be a list")

        if not len(pt_bounds) == 3:
            error("'pt_bounds' must have exactly 3 elements")
        else:
            for item in pt_bounds:
                if not is_int(item):
                    error("'pt_bounds' element not an integer")

        self.pt_bounds = pt_bounds

    def set_mass_spectrum(self, andi_data):

        """
        @summary: Sets peak mass spectrum

        @param andi_data: An ANDI data object (such as IO.ANDI.ChemStation)
        @type andi_data: ChemStation
        @return: none
        @rtype: NoneType
        """
        
        if self.pt_bounds == None:
            error("pt_bounds not set for this peak")
        else:
            pt_apex = self.pt_bounds[1]

        self.mass_spectrum = andi_data.get_mass_spectrum_at_index(pt_apex)
        self.mass_list = andi_data.get_mass_list()

        # These two must be consistent, this is checked in IO.ANDI
        # upon reading the ANDI data file.  Check again:
        if not (len(self.mass_spectrum) == len(self.mass_list)):
            error("peak mass spectrum data inconsistent")

    def crop_mass_spectrum(self, mass_min, mass_max):

        """
        @summary: Crops mass spectrum

        @param mass_min: Minimum mass value
        @type mass_min: IntType or FloatType 
        @param mass_max: Maximum mass value
        @type mass_max: IntType or FloatType
        @return: none
        @rtype: NoneType
        """
        
        if self.mass_spectrum == None:
            error("mass spectrum not set for this peak")
        else:
            if (not is_number(mass_min)) or (not is_number(mass_max)):
                error("'mass_min' and 'mass_max' must be numbers")

        new_mass_list = []
        new_mass_spectrum = [] 

        for ii in range(len(self.mass_list)):
            mass = self.mass_list[ii]
            intensity =  self.mass_spectrum[ii]
            if (mass >= mass_min) and (mass <= mass_max):
                new_mass_list.append(mass)
                new_mass_spectrum.append(intensity) 

        if len(new_mass_list) < 10:
            print " WARNING: peak mass spectrum contains < 10 points"
 
        self.mass_list = new_mass_list
        self.mass_spectrum = new_mass_spectrum

        if len(new_mass_list) == 0:
            error("mass spectrum empty")
        else:
            print " [ Mass spectrum cropped to %d points:" % \
                    (len(new_mass_list)),
            print "masses %d to %d ]" % (min(new_mass_list), 
                    max(new_mass_list))

    def set_peak_tag(self, tag):

        """
        @summary: Sets the peak tag

        @param tag: The supplied peak tag is either RF- type
            tag (such as "RF-SI") or "BLANK"
        @type tag: StringType
        @return: none
        @rtype: NoneType
        """

        peak_tag = None

        if tag != None:
            if tag[:3] == "rf-" and tag[3:] != "":
                peak_tag = tag 
            elif tag == "blank":
                peak_tag = tag
            else:
                error("incorrect reference peak tag '%s'" % (tag))

        self.tag = peak_tag

