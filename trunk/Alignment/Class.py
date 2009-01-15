"""
Classes for full matrix alignment by dynamic programming
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

import Function

import numpy

class Alignment(object):

    """
    @summary: Models an alignment of peak lists

    @param object: An experiment list
    @type object: pyms.Experiment.Class.Experiment

    @author: Qiao Wang
    @author: Vladimir Likic
    """

    def __init__(self, expr):

        """
        @summary: Initialize the alignment object

        @param expr: The experiment to be converted into an alignment object
        @type expr: pyms.Experiment.Class.Experiment

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        if expr == None:
            self.scanpos = []
            self.scangap = []
            self.similarity = None
        else:
            self.scanpos = []
            self.scangap = []
            data = expr.get_intensity_matrix()
            for i in range(len(data)):
                self.scanpos.append([])
                self.scanpos[i].append((numpy.array(data[i], dtype='d')).tolist())
                self.scangap.append([])
                self.scangap[i].append(1)
            self.similarity = None

    def __len__(self):

        """
        @summary: Returns the length of the alignment, defined as the number of
            peak positions in the alignment

        @return: the length of the alignment
        @rtype: IntType

        @author: Qiao Wang
        @author: Vladimir Likic
        """

        return len(self.scanpos)

class PairwiseAlignment(object):

    """
    @summary: Models pairwise alignment of alignments

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    print "Pairwise alignment started"

    def __init__(self, data1, data2, gap):

        """
        @param data1: First Alignment
        @type data1: pyms.Alignment.Class.Alignment
        @param data2: Second Alignment
        @type data2: pyms.Alignment.Class.Alignment
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        self.algt1 = data1
        self.algt2 = data2
        self.gap = gap
        self.sim_matrix = self._sim_matrix(self.algt1, self.algt2, self.gap)

    def _sim_matrix(self, algt1, algt2, gap):

        """
        @summary: Calculates the similarity matrix for the set of alignments

        @param algts: A list of alignments
        @type algts: ListType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @return: Similarity matrix
        @rtype: numpy.ndarray

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        sim_matrix = numpy.zeros((2, 2), dtype='f')
        ma = Function.align(algt1, algt2, gap)
        sim_matrix[0, 1] = sim_matrix[1, 0] = ma.similarity

        return sim_matrix