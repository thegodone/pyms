"""Class.py
 Module Class in metab.Peak
 Contains Peak-related classes.
"""

 #############################################################################
 #                                                                           #
 #    META-B software for processing of metabolomic mass-spectrometry data   #
 #    Copyright (C) 2005-6 Vladimir Likic                                    #
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

import sys, types, math
import numpy

from metab.Utils.Error import error
from metab.Utils.Utils import *

class Peak:

    """class Peak

    A class representing a signal peak. A peak object is
    initialised with retention time and raw peak area.

    def __init__(self, rt, raw_area)

    @param rt A number. Retention time.
    @param raw_area A number. Raw peak area.
    """

    def __init__(self, rt, raw_area):

        # Retention time and raw_area are mandatory
        self.rt = float(rt)
        self.raw_area = float(raw_area)

        # Other attributes
        self.intensity = None 
        self.pt_bounds = None
        self.rt_bounds = None
        self.norm_area = None 
        self.confidence = None
        self.mass_spectrum = None 
        self.mass_list = None 
        self.tag = None

    def set_intensity(self, intensity):

        """set_intensity(intensity)

        Sets peak intensity.

        @param intensity A number
        """

        if not is_number(intensity):
            error("'intensity' must be a number")

        self.intensity = float(intensity)

    def set_pt_bounds(self, pt_bounds):

        """set_pt_bounds(pt_bounds)

        Sets peak boundaries in points.

        @param pt_bounds A list. Contains left, apex, and right
            peak boundary (points).
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

    def set_rt_bounds(self, rt_bounds):

        """set_rt_bounds(rt_bounds)

        Sets peak retention time boundaries.

        @param rt_bounds A list. Contains left, apex, and right peak
            boundary (retention time).
        """

        if not is_list(rt_bounds):
            error("'rt_bounds' must be a list")

        if not len(rt_bounds) == 3:
            error("'rt_bounds' must have exactly 3 elements")
        else:
            for item in rt_bounds:
                if not is_number(item):
                    error("'rt_bounds' element not a number")

        if not math.fabs(rt_bounds[1] - self.rt) < 1e-7:
            error("rt_bounds[1] does not match peak retention time")

        self.rt_bounds = rt_bounds

    def set_norm_area(self, norm_area):

        """set_norm_area(self, norm_area)

        Sets peak normalised area. 

        @param norm_area. A number
        """

        if not is_number(norm_area):
            error("'norm_area' must be a number")

        self.norm_area = float(norm_area)

    def set_confidence(self, confidence):

        """set_confidence(self, confidence)

        Sets peak confidence level.

        @param confidence. A number between 0 and 1.
        """

        if not is_number(confidence):
            error("peak confidence must be a number")
        else:
            if confidence < 0.0 or confidence > 1.0:
                error("peak confidence must be between 0 and 1")

        self.confidence = float(confidence)

    def set_mass_spectrum(self, andi_data):

        """set_mass_spectrum(self, andi_data)

        Sets peak mass spectrum.

        @param andi_data An ANDI data object (such as IO.ANDI.ChemStation)
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

        """crop_mass_spectrum(self, mass_min, mass_max)
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

        """set_peak_tag(tag)

        Sets the peak tag.

        @param A string. The supplied peak tag is either RF- type
            tag (such as "RF-SI") or "BLANK".
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

