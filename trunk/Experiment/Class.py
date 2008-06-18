"""Class.py
 Module Class in pyms.Experiment
 Provides class Experiment.
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
from pyms.Utils.Utils import * 

from pyms import Peak

def cmp_peak_area(peak1, peak2):

    """cmp_peak_area(peak1, peak2)

    Compares peak normalised area.

    @param peak1 A peak object.
    @param peak2 A peak object.
    """

    return cmp(peak1.norm_area, peak2.norm_area)

class Experiment:

    """class Experiment

    A class representing a single experiment. 
    """

    def __init__(self, expr_code, peaks):

        """__init__(self, expr_code, peaks)

        @param expr_code A string. Code-name for the experiment
        @param peaks A list. A list of peak objects. 
        """

        self.expr_code = expr_code
        self.peaks = peaks
        self.ref_peak = None

    def set_ref_peak(self, ref_peak_tag):

        """set_ref_peak(self, ref_peak_tag)
        Finds the reference peak, sets the reference peak, and removes
        it from the list of peaks.

        @param ref_peak_tag A string. Tag string for the chosen reference
            peak. For example, two potential reference peaks may be
            annotated in the peak file, "RT-SI" and "RT-NV". If the first
            is to be used as the reference peak, ref_peak_tag should be
            set to "si".
        """

        ref_peaks = {}

        for peak in self.peaks:
            if (peak.tag != None) and (peak.tag[:3] == "rf-"):
                ref_tag = peak.tag[3:]
                if ref_tag in ref_peaks.keys():
                    error("multiple reference peaks with the tag '%s'" \
                    % (peak.tag))
                else:
                    ref_peaks[ref_tag] = peak

        if len(ref_peaks.keys()) == 0:
            print " Experiment: %s :" % (self.expr_code),
            error("no reference peaks found")

        ref_peak = None

        for ref_tag in ref_peaks.keys():
             if ref_tag == ref_peak_tag:
                 ref_peak = ref_peaks[ref_tag]

        if ref_peak == None:
            print " Experiment: %s :" % (self.expr_code),
            error("specified reference peak '%s' not found" % \
                    (ref_peak_tag))
        else:
            self.ref_peak = ref_peak
            print "\t[ Reference peak found: '%s' @ %.3f s" % \
                    (self.ref_peak.tag, self.ref_peak.rt)

        # remove all reference peaks from the peak list
        for ref_tag in ref_peaks.keys():
            rm_ref_peak = ref_peaks[ref_tag]
            print "\t  [ Removing reference peak '%s' @ %.3f s ]" % \
                (rm_ref_peak.tag, rm_ref_peak.rt)
            self.peaks.remove(rm_ref_peak)

    def remove_blank_peaks(self, peak_remove_tag):

        """remove_blank_peaks(self, peak_remove_tag)
        Removes peaks which whose tag is 'peak_remove_tag'.

        @param peak_remove_tag A string. Peaks tagged by this tag will
            be discarded. Normally, peak_remove_tag equals "blank".
        """

        remove_list = []
        for peak in self.peaks:
            if peak.tag == peak_remove_tag:
                remove_list.append(peak)

        for peak in remove_list:
            self.peaks.remove(peak)
            print "\t[ Designated blank peak at %.3f s removed ]" \
                        % (peak.rt)

        if len(remove_list) == 0:
            print "\t[no peaks to be removed]"

    def normalise_peaks(self, to_reference=True):

        """normalise_peaks(to_reference)
        Normalizes peak areas.
        """

        if to_reference:
            ref_area = self.ref_peak.raw_area
        else:
            ref_area = 1.0

        for peak in self.peaks:
            peak.norm_area = float(peak.raw_area)/float(ref_area)

    def sele_top_peaks(self, n):

        """sele_top_peaks(self, n)

        Selects n strongest peaks by normalised area, discards others.

        @param n integer.
        """

        self.peaks.sort(cmp_peak_area)
        self.peaks.reverse()
        self.peaks = self.peaks[:n]

        print " -> %d strongest peaks selected" % ( len(self.peaks) )

    def purge_peaks(self, norm_area_threshold):

        """purge_peaks(self, norm_area_threshold)

        Purge peaks which are below the threshold expressed as the
        fraction of the _normalized_ reference peak area (i.e. this
        is exactly a normalized peak area peak.norm_area).

        @param norm_area_threshold A float.
        """

        purge_list = []
        for peak in self.peaks:
            if peak.norm_area < norm_area_threshold:
                purge_list.append(peak)

        for peak in purge_list:
            self.peaks.remove(peak)

        print " Experiment %s: %d peaks purged (below threshold=%.2f)" % \
                (self.expr_code, len(purge_list), norm_area_threshold )

    def scale_peaks(self, scale_factor):

        """scale_peaks(self, scale_factor)

        Scales all normalised peak areas by a factor. Typically used to
        adjust DNA normalisation scaling. 

        @param scale_factor A number.
        """

        if not is_number(scale_factor):
            error("scale factor not a number")

        f = float(scale_factor)

        print " -> Scaling peak normalised areas with %.2f" % ( f )

        for peak in self.peaks:
            peak.norm_area = peak.norm_area*f

    def sele_rt_range(self, rt_range):

        """sele_rt_range(self, rt_range)

        Discards all peaks with retention times outside the specified
        range.

        @param rt_range A list. Contains two numbers [rt_min, rt_max].
        """

        print " -> Selection by retention time (from %s to %s):" % \
                (rt_range[0], rt_range[1]),

        peaks_sele = Peak.List.Utils.sele_peaks_by_rt(self.peaks, rt_range)
        self.peaks = peaks_sele


