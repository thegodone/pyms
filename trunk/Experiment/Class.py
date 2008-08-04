"""
Models an GC-MS or LC-MS experiment represented by a list of signal peaks
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
from pyms import Peak
from pyms.Peak.List.Utils import sele_peaks_by_rt 

class Experiment:

    """
    @summary: Models an experiment object

    @author: Vladimir Likic
    """

    def __init__(self, expr_code, peaks):

        """
        @param expr_code: Code-name for the experiment
        @type expr_code: StringType
        @param peaks: A list of peak objects
        @type: ListType
        """

        self.expr_code = expr_code
        self.peaks = peaks
        self.ref_peak = None

    def set_ref_peak(self, ref_peak_tag):

        """
        @summary: Sets the reference peak and removes it from the list of peaks

        Reference peaks are annotated with the RF-type tags. For example, two
        potential reference peaks may be annotated in the peak file, "RT-SI"
        and "RT-NV". To set the first peak as the reference peak, 'ref_peak_tag'
        should be set to "si".

        @param ref_peak_tag: Tag string for the chosen reference peak
        @type ref_peak_tag: StringType

        @return: none
        @rtype: NoneType
        """

        ref_peaks = {}

        for peak in self.peaks:
            if peak.tag != None and peak.tag[:3] == "rf-":
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
            print " [ Reference peak found: '%s' @ %.3f s ]" % \
                    (self.ref_peak.tag, self.ref_peak.rt)

        # remove all reference peaks from the peak list
        for ref_tag in ref_peaks.keys():
            rm_ref_peak = ref_peaks[ref_tag]
            print "  [ Removing reference peak '%s' @ %.3f s ]" % \
                (rm_ref_peak.tag, rm_ref_peak.rt)
            self.peaks.remove(rm_ref_peak)

    def remove_blank_peaks(self, peak_remove_tag="blank"):

        """
        @summary: Removes peaks whose tag equals 'peak_remove_tag'

        @param peak_remove_tag: Peaks tagged by this tag will be discarded
            from the peak list
        @type peak_remove_tag: StringType

        @return: none
        @rtype: NoneType
        """

        remove_list = []

        for peak in self.peaks:
            if peak.tag == peak_remove_tag:
                remove_list.append(peak)

        for peak in remove_list:
            self.peaks.remove(peak)
            print "\t[ Designated blank peak at %.3f s removed ]" \
                        % ( peak.rt )

        if len(remove_list) == 0:
            print "\t[no peaks to be removed]"

    def raw2norm_area(self, to_reference=False):

        """
        @summary: Normalises peak areas

        This method populates peak attributes 'norm_area'. If to_reference=True,
        peak areas are normalized by the raw area of the reference peak (assumed
        to be previously set). If to_reference=False, 'norm_area' is the same as
        'raw_area'.

        @param to_reference: True for normalization to reference peak, False
        otherwise
        @type to_reference: BooleanType

        @return: None
        @rtype: NoneType
        """

        if to_reference:
            ref_area = self.ref_peak.raw_area
        else:
            ref_area = 1.0

        for peak in self.peaks:
            peak.norm_area = float(peak.raw_area)/float(ref_area)

    def purge_peaks(self, norm_area_threshold=0.0):

        """
        @summary: Purge peaks which are below the threshold expressed as the
        fraction peak.norm_area

        @param norm_area_threshold: Threshold value
        @type norm_area_threshold: IntType or FloatType

        @return: None
        @rtype: NoneType
        """

        purge_list = []

        for peak in self.peaks:
            if peak.norm_area < norm_area_threshold:
                purge_list.append(peak)

        for peak in purge_list:
            self.peaks.remove(peak)

        print " Experiment %s: %d peaks purged (below threshold=%.2f)" % \
                (self.expr_code, len(purge_list), norm_area_threshold )
 
    def sele_rt_range(self, rt_range):

        """
        @summary: Discards all peaks which have the retention time outside
        the specified range

        @param rt_range: Min, max retention time given as a list [rt_min, rt_max]
        @type rt_range: ListType

        @return: none
        @rtype: NoneType
        """

        print " -> Selecting peaks by retention time (from %s to %s):" % \
                (rt_range[0], rt_range[1]),

        peaks_sele = sele_peaks_by_rt(self.peaks, rt_range)
        self.peaks = peaks_sele

