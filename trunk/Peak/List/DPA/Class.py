"""Class.py
"""

import numpy
import Pycluster

from pyms.Utils.Error import error 
from pyms.Utils.DP import dp 

from pyms import Experiment

import Function

class Alignment:

    """ Holds a list of Alignments """

    def __init__(self, experiment=None):

        """ Initialization

        @param experiment Provide this if we're making a singleton alignment
            from an Experiment (optional)
        """

        if experiment and isinstance(experiment, Experiment.Class.Experiment):
            self.alignments = [experiment.peaks[:]] # work on copy of the peaklist
            self.expr_codes = [experiment.expr_code]
            self.similarity = 0.

    def __len__(self):
        return len(self.alignments)

    def write_csv(self, rt_filename, area_filename, minutes=True):

        """ Writes the alignments into a CSV file

        @param rt_filename File to save retention times to
        @param area_filename File to save normalized areas to
        @param minutes True/False to save retention time in minutes
        """

        try:
            fp_rt = open(rt_filename, 'w')
            fp_area = open(area_filename, 'w')
        except IOError:
            error("Cannot open output file for writing")

        # write experiment headers
        header = '"' + '","'.join(self.expr_codes) + "\"\n"

        fp_rt.write(header)
        fp_area.write(header)

        # for each peak, write each alignment's peak and area value
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

    def transpose(self):

        """ Transposes the alignment lists (i.e. makes one row per peak match VS
            one row per experiment. For convenient processing by some routines. 
        """

        self.alignments = [ [r[col] for r in self.alignments]
                            for col in range(len(self.alignments[0]))
                          ]

    def filter_min_peaks(self, min_peaks):

        """ Filters away peak matches that have less than min_peaks

        @param min_peaks Minimum # of peaks
        """

        self.transpose()
        self.alignments =  [ x for x in self.alignments
                             if len(filter(None, x)) >= min_peaks
                           ]
        self.transpose()

class Tree:

    """ Holds a guide tree to be used with align_with_tree() """

    def __init__(self, reps, D, gap):
        
        """ Initialization

        @param reps A list of replicates. Can consist of Alignments or Experiements.
        @param D D parameter
        @param gap Gap parameter
        """

        self.D = D
        self.gap = gap

        # Handle trivial cases of 0 and 1
        if len(reps) == 0:
            return

        if len(reps) == 1:
            return

        n = len(reps)
        total = n * (n - 1) / 2

        # Step 1 - do pairwise alignments and make similarity matrix
        print " Calculating similarity matrix for %d replicates (D=%.2f, gap=%.2f)" % (n, D, gap)
        sim_matrix = numpy.zeros((n,n), dtype='f')
        for i in range(n - 1):
            for j in range(i + 1, n):
                alignment = Function.align(reps[i], reps[j], D, gap)
                sim_matrix[i,j] = sim_matrix[j,i] = alignment.similarity

                total = total - 1
                print " -> %d replicate pairs remaining" % total

	# Step 2 - change similarity matrix entries (i,j) to: max {matrix} - (i,j)
	sim_max = numpy.max(numpy.ravel(sim_matrix))
	sim_matrix = sim_max - sim_matrix
	for i in range(n):
		sim_matrix[i,i] = 0

        self.sim_matrix = sim_matrix

        # Step 3 - perform hierchical clustering
        self.tree = Pycluster.treecluster(distancematrix=sim_matrix, method='a')

