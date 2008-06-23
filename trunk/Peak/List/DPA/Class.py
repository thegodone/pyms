"""
Classes for alignment of peak lists by dynamic programming 
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
import Pycluster

from pyms.Utils.Error import error 
from pyms.Utils.DP import dp 

from pyms.Experiment.Class import Experiment

import Function

class Alignment:

    """
    @summary: Models an alignment of peak lists

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    def __init__(self, experiment=None):

        """
        @param experiment: Provide this to make a singleton alignment object from
            an experiment
        @type experiment: pyms.Experiment.Class.Experiment
        """

        if experiment and isinstance(experiment, Experiment):
            self.alignments = [experiment.peaks[:]] # work on copy of the peaklist
            self.expr_codes = [experiment.expr_code]
            self.similarity = 0.

    def __len__(self):

        """
        @summary: Returns the length of the alignment
        """

        return len(self.alignments)

    def transpose(self):

        """
        @summary: Transposes the alignment

        @author: Woon Wai Keen
        """

        self.alignments = [ [r[col] for r in self.alignments]
                            for col in range(len(self.alignments[0]))
                          ]

    def filter_min_peaks(self, min_peaks):

        """
        @summary: Filters alignment positions that have less than 'min_peaks'

        @param min_peaks: Minimum number of peaks required for the alignment
            position to survive filtering
        @type min_peaks: IntType

        @author: Woon Wai Keen
        """

        self.transpose()

        self.alignments =  [ x for x in self.alignments
                             if len(filter(None, x)) >= min_peaks
                           ]
        self.transpose()

    def write_csv(self, rt_file_name, area_file_name, minutes=True):

        """
        @summary: Writes the alignment to CSV files

        This function writes two files: one containing the alignment of peak
        retention times and the other containing the alignment of peak areas.

        @param rt_file_name: The name for the retention time alignment file
        @type rt_file_name: StringType
        @param area_file_name: The name for the areas alignment file
        @type area_file_name: StringType
        @param minutes: An indicator whether to save retention times in
            minutes. If False, retention time will be saved in seconds 
        @type minutes: BooleanType
        @rtype: NoneType

        @author: Woon Wai Keen
        """

        try:
            fp_rt = open(rt_file_name, 'w')
            fp_area = open(area_file_name, 'w')
        except IOError:
            error("Cannot open output file for writing")

        # write experiment headers
        header = '"' + '","'.join(self.expr_codes) + "\"\n"

        fp_rt.write(header)
        fp_area.write(header)

        # for each alignment position write each alignment's peak and area
        for peak_idx in range(len(self.alignments[0])):
            rts = []
            areas = []
            for align_idx in range(len(self.alignments)):
                peak = self.alignments[align_idx][peak_idx]
                if peak is not None:
                    if minutes:
                        rtmin = peak.rt/60.0
                    else: 
                        rtmin = peak.rt
                    rts.append('%.3f' % rtmin)
                    areas.append('%.4f' % peak.norm_area)
                else:
                    rts.append('NA')
                    areas.append('NA')

            fp_rt.write(",".join(rts) + "\n")
            fp_area.write(",".join(areas) + "\n")

        fp_rt.close()
        fp_area.close()

class Tree:

    """
    @summary: Models pairwise alignment on the set of experiments, and builds
    a guide tree

    This class expects a set of experiments (ie. the 'Experiment' objects) given
    as a list, and the alignment parameters.

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    def __init__(self, reps, D, gap):
        
        """
        @param reps: A set of experiments 
        @type reps: A list of Experiment instances
        @param D: Retention time tolerance parameter for pairwise alignments
        @type D: FloatType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        self.D = D
        self.gap = gap

        self.sim_matrix = self._sim_matrix(self.D, self.gap, reps)

        if self.sim_matrix == None:
            self.dist_matrix = None
            self.tree = None
        else:
            self.dist_matrix = self._dist_matrix(self.sim_matrix)
            self.tree = self._guide_tree(self.dist_matrix)

    def _sim_matrix(self, D, gap, reps):

        """
        @summary: calculates the similarity matrix for the set of replicate
        experiments

        @param D: Retention time tolerance parameter for pairwise alignments
        @type D: FloatType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType
        @param reps: A list of experiments 
        @type reps: A list of Experiment instances

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        # handle trivial cases
        if len(reps) in [0, 1]: return None

        n = len(reps)
        total_n = n * (n - 1) / 2

        print " Calculating pairwise alignments for %d replicates (D=%.2f, gap=%.2f)" % \
                (n, D, gap)

        sim_matrix = numpy.zeros((n,n), dtype='f')

        for i in range(n - 1):
            for j in range(i + 1, n):
                alignment = Function.align(reps[i], reps[j], D, gap)
                sim_matrix[i,j] = sim_matrix[j,i] = alignment.similarity
                total_n = total_n - 1
                print " -> %d replicate pairs remaining" % total_n

        return sim_matrix

    def _dist_matrix(self, sim_matrix):

        """
        @summary: Converts similarity matrix into a distance matrix

        @param sim_matrix: the similarity matrix
        @type sim_matrix: numpy.ndarray

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        # change similarity matrix entries (i,j) to max{matrix}-(i,j)
        sim_max = numpy.max(numpy.ravel(sim_matrix))
        dist_matrix = sim_max - sim_matrix

        # set diagonal elements of the similarity matrix to zero
        for i in range(len(dist_matrix)):
            dist_matrix[i,i] = 0

        return dist_matrix

    def _guide_tree(self, dist_matrix):

        """
        @summary: Build a guide tree from the distance matrix

        @param dist_matrix: distance matrix
        @type sim_matrix: numpy.ndarray

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        n = len(dist_matrix)

        print " Clustering %d pairwise alignments." % (n*(n-1)),

        tree = Pycluster.treecluster(distancematrix=dist_matrix, method='a')

        print "Done"

        return tree


