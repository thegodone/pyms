"""
Functions for full matrix alignment by dynamic programming
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

import Class

import numpy

from Class import PairwiseAlignment
from pyms.Utils.DP import dp
from pyms.Utils.Utils import is_list

def fma(data1, data2, Gw):

    """
    @summary: Align 2 data sets with full matrix

    @param data1: First input of GC-MS data
    @type data1: pyms.IO.ANDI.Class.ChemStation
    @param data2: Second input of GC-MS data
    @type data2: pyms.IO.ANDI.Class.ChemStation

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    print "ANDI-MS data 1 filename:", data1.get_filename()
    print "ANDI-MS data 2 filename:", data2.get_filename()
    # print the name of the ANDI-MS file
    im1 = data1.get_intensity_matrix()
    im2 = data2.get_intensity_matrix()
    # get the entire intensity matrix
    print "Dimensions of the intensity matrix 1 are:", len(im1), "x", len(im1[0])
    print "Dimensions of the intensity matrix 2 are:", len(im2), "x", len(im2[0])
    E1 = [im1, im2]
    F1 = exprl2alignment(E1)
    PairwiseAlignment(F1, Gw)

def exprl2alignment(exprl):

    """
    @summary: Converts the experiments into alignment objects

    @param exprl: An experiments list
    @type exprl: ListType

    @return: A list of alignments objects
    @rtype: ListType

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    if not is_list(exprl):
        error("the argument is not a list")

    algts = []
    for item in exprl:
        algt = Class.Alignment(item)
        algts.append(algt)

    print "Alignment objects created"

    return algts

def align(a1, a2, gap):

    """
    @summary: Aligns two alignments

    @param a1: The first alignment
    @type a1: pyms.Alignment.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Alignment.Class.Alignment
    @param gap: Gap penalty
    @type gap: FloatType

    @return: Aligned alignments
    @rtype: pyms.Alignment.Class.Alignment

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    # calculate score matrix for two alignments
    M = score_matrix(a1, a2)
    # run dynamic programming
    result = dp(M, gap)
    # make composite alignment from the results
    ma = merge_alignments(a1, a2, result['trace'])
    # calculate the similarity score
    ma.similarity = alignment_similarity(result['trace'], M, gap)
    return ma

def score_matrix(a1, a2):

    """
    @summary: Calculates the score matrix between two alignments

    @param a1: The first alignment
    @type a1: pyms.Alignment.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Alignment.Class.Alignment

    @return: A substitution matrix for dynamic programming
    @rtype: numpy.ndarray

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    score_matrix = numpy.zeros((len(a1.scanpos), len(a2.scanpos)))
    print "Initializing the score matrix of size:", score_matrix.shape
    print "Processing the score matrix..."

    row = 0
    col = 0

    for pos1 in a1.scanpos:
        for pos2 in a2.scanpos:
            score_matrix[row][col] = position_similarity(pos1, pos2)
            col = col + 1
        row = row + 1
        col = 0

    print "Score matrix finished"

    return score_matrix

def position_similarity(pos1, pos2):

    """
    @summary: Calculates the similarity between each alignment pairs of positions

    @param pos1: The position of the first alignment
    @type pos1: ListType
    @param pos2: The position of the second alignment
    @type pos2: ListType

    @return: The similarity value for the current position between each pair
    @rtype: FloatType

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    score = 0.0
    mass_spect1 = numpy.array(pos1[0], dtype='d')
    mass_spect2 = numpy.array(pos2[0], dtype='d')
    mass_spect1_sum = numpy.sum(mass_spect1 ** 2, axis=0)
    mass_spect2_sum = numpy.sum(mass_spect2 ** 2, axis=0)
    top = numpy.dot(mass_spect1, mass_spect2)
    bot = numpy.sqrt(mass_spect1_sum * mass_spect2_sum)
    cos = top / bot
    score = score + (1.0 - cos)

    return score

def merge_alignments(A1, A2, traces):

    """
    @summary: Merges two alignments with gaps added in from DP traceback

    @param A1: First alignment
    @type A1: pyms.Alignment.Class.Alignment
    @param A2: Second alignment
    @type A2: pyms.Alignment.Class.Alignment
    @param traces: DP traceback
    @type traces: ListType

    @return: A single alignment from A1 and A2
    @rtype: pyms.Alignment.Class.Alignment

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    # Create object to hold new merged alignment and fill in its expr_codes
    print "Merging the alignments started"
    ma = Class.Alignment(None)

    # create empty lists of dimension |A1| + |A2|
    idx1 = idx2 = 0

    # trace can either be 0, 1, or 2
    # if it is 0, there are no gaps. otherwise, if it is 1 or 2,
    # there is a gap in A2 or A1 respectively.
    count = 0
    for trace in traces:

        if trace == 0:
            ma.scanpos.append([])
            ma.scangap.append([])
            for i in range(len(A1.scanpos[idx1])):
                ma.scanpos[count].append(A1.scanpos[idx1][i])
                ma.scangap[count].append(1)
            for j in range(len(A2.scanpos[idx2])):
                ma.scanpos[count].append(A2.scanpos[idx2][j])
                ma.scangap[count].append(1)
            idx1 = idx1 + 1
            idx2 = idx2 + 1

        elif trace == 1:
            ma.scanpos.append([])
            ma.scangap.append([])
            for i in range(len(A1.scanpos[idx1])):
                ma.scanpos[count].append(A1.scanpos[idx1][i])
                ma.scangap[count].append(1)
            for j in range(len(A2.scanpos[idx2])):
                ma.scanpos[count].append([])
                ma.scangap[count].append(0)
            idx1 = idx1 + 1

        elif trace == 2:
            ma.scanpos.append([])
            ma.scangap.append([])
            for i in range(len(A1.scanpos[idx1])):
                ma.scanpos[count].append([])
                ma.scangap[count].append(0)
            for j in range(len(A2.scanpos[idx2])):
                ma.scanpos[count].append(A2.scanpos[idx2][j])
                ma.scangap[count].append(1)
            idx2 = idx2 + 1

        count = count + 1

    print "Start processing artificial time points"
    fp = open('Artificial_value.txt', 'w')
    for i in range(len(ma.scanpos[0])):
        for j in range(len(ma.scanpos)):
            if len(ma.scanpos[j][i]) == 0:
                start = j
                fp.write("-----------------------------------------------\n")
                fp.write("Gaps founded, Start at: " + str(start) + " , ")
                end = start
                while len(ma.scanpos[end + 1][i]) == 0:
                    end = end + 1
                fp.write("End at: " + str(end) + " .\n")
                fp.write("The details of artificial values of masspectrum:\n")
                if j == 0:
                    for l in range(len(ma.scanpos[0][i])):
                        slope = (ma.scanpos[end + 1][i][l]-0) / (end-start + 1)
                        ma.scanpos[0][i].append(0.0)
                        for k in range(end-start):
                            ma.scanpos[k + 1][i].append(0.0 + slope * (k + 1))
                elif end == len(ma.scanpos)-1:
                    for l in range(len(ma.scanpos[0][i])):
                        slope = (ma.scanpos[start-1][i][l]-0.0) / (end-start + 1)
                        for k in range(end-start + 1):
                            ma.scanpos[start + k].append(ma.scanpos[start-1] - slope * (k + 1))
                else:
                    for l in range(len(ma.scanpos[0][i])):
                        slope = (ma.scanpos[end + 1][i][l]-ma.scanpos[start-1][i][l]) / (end-start + 2)
                        fp.write("Gap point:  Slope: " + str(slope))
                        fp.write(" Left: " + str(ma.scanpos[start-1][i][l]))
                        for k in range(end-start + 1):
                            ma.scanpos[start + k][i].append(ma.scanpos[start-1][i][l] + slope * (k + 1))
                            fp.write(" Inter: " + str(ma.scanpos[start + k][i][l]))
                        fp.write(" Right: " + str(ma.scanpos[end + 1][i][l]) + "\n")
                    fp.write("-----------------------------------------------\n")
    fp.close()

    print "Artificial time points filled with interpolated mass spectrum value"
    print "Generating the final 2-alignment..."
    fp = open('Final_2-alignment.txt', 'w')
    for i in range(len(ma.scanpos)):
        for j in range(len(ma.scanpos[0])):
            a = 0
            for k in range(len(ma.scanpos[i][j])):
                a = a + ma.scanpos[i][j][k]
            fp.write(str(a) + " ")
        fp.write("\n")
    print "Done"
    fp.close()
    print "Generating the gap indicator matrix..."
    fp = open('Gaps_indicator.txt', 'w')
    for i in range(len(ma.scangap)):
        for j in range(len(ma.scangap[0])):
            fp.write(str(ma.scangap[i][j]) + " ")
        fp.write("\n")
    print "Done"
    fp.close()
    print "Merging complete"
    return ma

def alignment_similarity(traces, score_matrix, gap):

    """
    @summary: Calculates similarity score between two alignments

    @param traces: Traceback from DP algorithm
    @type traces: ListType
    @param score_matrix: Score matrix of the two alignments
    @type score_matrix: numpy.ndarray
    @param gap: Gap penalty
    @type gap: FloatType

    @return: Similarity score (i.e. more similar => higher score)
    @rtype: FloatType
    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    score_matrix = 1.0 - score_matrix
    similarity = 0.0
    idx1 = idx2 = 0
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