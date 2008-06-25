"""
Functions for peak alignment by dynamic programming
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

import copy

import numpy

from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import is_list 
from pyms.Utils.DP import dp

from pyms.Experiment.Class import Experiment

import Class
import Utils

def align_with_tree(T, min_peaks=1):

    """
    @summary: Aligns a list of alignments using the supplied guide tree

    @param T: The pairwise alignment object
    @type: pyms.Peak.List.DPA.Class.PairwiseAlignment 
    @return: The final alignment consisting of aligned input alignments
    @rtype: pyms.Peak.List.DPA.Class.Alignment
    """

    print " Aligning %d items with guide tree (D=%.2f, gap=%.2f)" % \
            (len(T.algts), T.D, T.gap)

    # For everything else, we align according to the guide tree provided by
    # Pycluster. From Pycluster documentation:
    #   Each item and subnode is represented by an integer. For hierarchical
    #   clustering of n items, we number the original items {0, ... , n-1},
    #   nodes are numbered {-1, ... , -(n-1)}. Note that the number of nodes
    #   is one less than the number of items.

    # extend As to length 2n to hold the n items, n-1 nodes, and 1 root
    As = copy.deepcopy(T.algts) + [ None for _ in range(len(T.algts)) ]

    # align the alignments into positions -1, ... ,-(n-1)
    total = len(T.tree)
    index = 0

    for node in T.tree:
        index = index - 1
        As[index] = align(As[node.left], As[node.right], T.D, T.gap)
        total = total - 1
        print " -> %d item(s) remaining" % total

    # the final alignment is in the root. Filter min peaks and return
    As[index].filter_min_peaks(min_peaks)

    return As[index]

def expr2alignment(exprl):

    """
    @summary: Converts experiment(s) into alignment(s)

    @param exprl: Experiment(s) to be converted into an alignment object(s)
    @type exper: pyms.Experiment.Class.Experiment instance or list of such instances

    @author: Vladimir Likic
    """

    if isinstance(exprl, Experiment):

        algt = Class.Alignment(exprl)

        return algt

    elif is_list(exprl):

        algts = []

        for item in exprl:

            if not isinstance(item, Experiment):
                error("list items must be 'Experiment' instances")
            else:
                algt = Class.Alignment(item)

            algts.append(algt)

        return algts

    else:
        error("the argument must be 'Experiment' instance, or list of such instances")

def align(a1, a2, D, gap):

    """ 
    @summary: Aligns two alignments

    @param a1: The first alignment
    @type a1: pyms.Peak.List.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Peak.List.Class.Alignment
    @param D: Retention time tolerance
    @type D: FloatType
    @param gap: Gap penalty
    @type D: FloatType

    @return: Aligned alignments
    @rtype: pyms.Peak.List.Class.Alignment

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    # calculate score matrix for two alignments
    M = score_matrix(a1, a2, D)

    # run dynamic programming
    result = dp(M, gap)

    # make composite alignment from the results
    ma = merge_alignments(a1, a2, result['trace'])

    # calculate the similarity score
    ma.similarity = alignment_similarity(result['trace'], M, gap)

    return ma

def merge_alignments(A1, A2, traces):

    """
    @summary: Merges two alignments with gaps added in from DP traceback

    @param A1 First alignment
    @param A2 Second alignment
    @param traces DP traceback

    @return A single alignment from A1 and A2
    """

    # Create object to hold new merged alignment and fill in its expr_codes
    ma = Class.Alignment(None)
    ma.expr_code = A1.expr_code + A2.expr_code

    # create empty lists of dimension |A1| + |A2|
    merged = [ [] for _ in range(len(A1.peakpos) + len(A2.peakpos)) ]

    A1 = A1.peakpos
    A2 = A2.peakpos

    idx1 = idx2 = 0

    # trace can either be 0, 1, or 2
    # if it is 0, there are no gaps. otherwise, if it is 1 or 2,
    # there is a gap in A2 or A1 respectively.

    for trace in traces:

        if trace == 0:

            for i in range(len(A1)):
                merged[i].append(A1[i][idx1])

            for j in range(len(A2)):
                merged[1+i+j].append(A2[j][idx2])

            idx1 = idx1 + 1
            idx2 = idx2 + 1

        elif trace == 1:

            for i in range(len(A1)):
                merged[i].append(A1[i][idx1])

            for j in range(len(A2)):
                merged[1+i+j].append(None)

            idx1 = idx1 + 1

        elif trace == 2:

            for i in range(len(A1)):
                merged[i].append(None)

            for j in range(len(A2)):
                merged[1+i+j].append(A2[j][idx2])

            idx2 = idx2 + 1

    ma.peakpos = merged

    # sort according to average peak
    ma.transpose()
    ma.peakpos.sort(Utils.alignment_compare)
    ma.transpose()

    return ma

def alignment_similarity(traces, score_matrix, gap):

    """ Calculates similarity score between two alignments (new method)

    @param traces Traceback from DP algorithm
    @param score_matrix Score matrix of the two alignments
    @param gap Gap penalty

    @return Similarity score (i.e. more similar => higher score)
    """

    score_matrix = 1. - score_matrix
    similarity = 0.
    idx1 = idx2 = 0

    # Trace can either be 0, 1, or 2
    # If it is 0, there is a match and we add to the sum the score between
	# these two aligned peaks. 
	# 
	# Otherwise, if it is 1 or 2, and there is a gap in A2 or A1
	# respectively. We then subtract the gap penalty from the sum.
    for trace in traces:
        if trace == 0:
	    similarity = similarity + score_matrix[idx1][idx2]
            idx1 = idx1 + 1
            idx2 = idx2 + 1
        elif trace == 1:
	    similarity = similarity - gap
            idx1 = idx1 + 1
        elif trace == 2:
	    similarity = similarity - gap
            idx2 = idx2 + 1

    return similarity

################################################################

def score_matrix(a1, a2, D):

    """
    @summary: Calculates the score matrix between two alignments

    @param a1: The first alignment
    @type a1: pyms.Peak.List.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Peak.List.Class.Alignment
    @param D: Retention time tolerance
    @type D: FloatType

    @return: Aligned alignments
    @rtype: pyms.Peak.List.Class.Alignment

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    score = numpy.zeros((len(a1), len(a2)))
    times = numpy.zeros((len(a1), len(a2)))

    # for all pairs of peak positions calculate the distance between them
    for peaks1 in a1.peakpos:
        for peaks2 in a2.peakpos:
            score_matrix, bin_matrix = score_matrix_pair(peaks1, peaks2, D)
            score = score + score_matrix
            times = times + bin_matrix

    # divide score_matrix by # of times the pair's score was calculated
    score = score/times

    return score

def score_matrix_pair(peaks1, peaks2, D):

    """
    @summary: Calculates the score matrix between two peaklists

    @param peaks1: The first list of peaks
    @param peaks2: The second list of peaks
    @param D: Retention time tolerance

    @return A tuple (score matrix, binary matrix)
    """

    # Initialize gap matrix, gap1 and gap2. This is the identity matrix with
    # zero rows/columns inserted where there's a gap in peaks1/peaks2
    f_peaks1 = filter(None, peaks1)
    f_peaks2 = filter(None, peaks2)

    gap1 = numpy.zeros([len(peaks1), len(f_peaks1)])
    gap2 = numpy.zeros([len(f_peaks2), len(peaks2)])

    col = 0
    for row in range(len(peaks1)):
        if peaks1[row] is not None:
            gap1[row][col] = 1
            col = col + 1

    row = 0
    for col in range(len(peaks2)):
        if peaks2[col] is not None:
            gap2[row][col] = 1
            row = row + 1

    # Create mass spectra vectors and retention time vectors
    ms1 = numpy.transpose(numpy.array([ p.mass_spectrum for p in f_peaks1 ], dtype='d'))
    ms2 = numpy.transpose(numpy.array([ p.mass_spectrum for p in f_peaks2 ], dtype='d'))

    rt1 = [ p.rt for p in f_peaks1 ]
    rt2 = [ p.rt for p in f_peaks2 ]

    # Do cosine angle transformation and modulation between all pairs of peaks
    cos_matrix = 1. - cos_ndp(ms1, ms2)
    mod_matrix = modulate(rt1, rt2, D)

    score_matrix = 1. - (cos_matrix * mod_matrix)

    # Insert gaps into matrix where there were Nones before, and also create
    # binary matrix, used to see which peak pairs had calculations on them
    score_matrix = numpy.dot(numpy.dot(gap1, score_matrix), gap2)
    bin_matrix = numpy.dot(numpy.dot(gap1, numpy.ones([len(f_peaks1), len(f_peaks2)])), gap2)

    return score_matrix, bin_matrix

def cos_ndp(ms1, ms2):

    """ Calculates cosine distance between two mass spectra

    @param ms1 First mass spectrum matrix
    @param ms2 Second mass spectrum matrix

    @return Cosine distance between ms1 and ms2
    """
    sum_ms1squared = numpy.sum(ms1 ** 2, axis=0)
    sum_ms2squared = numpy.sum(ms2 ** 2, axis=0)
 
    # perform matrix multiplication on above
    top = numpy.dot(numpy.transpose(ms1), ms2)

    sum_ms1squared = numpy.reshape(sum_ms1squared, [-1, 1])
    sum_ms2squared = numpy.reshape(sum_ms2squared, [1, -1])

    bot = numpy.sqrt(numpy.dot(sum_ms1squared, sum_ms2squared))

    cs = 1. - (top/bot)

    return cs

def modulate(rt1, rt2, D):

    """ Calculates gaussian distance between two peak retention time lists

    @param rt1 First retention time list
    @param rt2 Second retention time list
    @param D D

    @return Gaussian similarity matrix between rt1 and rt2
    """

    M = numpy.zeros(shape=(len(rt1), len(rt2)), dtype='d')

    for i in range(len(rt1)):
        for j in range(len(rt2)):
            M[i][j] = numpy.exp(-((rt1[i]-rt2[j]) / D)**2 / 2.)

    return M

