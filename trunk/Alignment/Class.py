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

import Pycluster

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
            for i in range(len(expr)):
                self.scanpos.append([])
                self.scanpos[i].append((numpy.array(expr[i], dtype='d')).tolist())
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

    def __init__(self, algts, gap):

        """
        @param algts: A list of alignments
        @type algts: ListType
        @param gap: Gap parameter for pairwise alignments
        @type gap: FloatType

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        self.algts = algts
        self.gap = gap
        self.sim_matrix = self._sim_matrix(algts, gap)
        self.dist_matrix = self._dist_matrix(self.sim_matrix)

    def _sim_matrix(self, algts, gap):

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
        ma = Function.align(algts[0], algts[1], gap)
        sim_matrix[0, 1] = sim_matrix[1, 0] = ma.similarity

        return sim_matrix

    def _dist_matrix(self, sim_matrix):

        """
        @summary: Converts similarity matrix into a distance matrix

        @param sim_matrix: The similarity matrix
        @type sim_matrix: numpy.ndarray
        @return: Distance matrix
        @rtype: numpy.ndarray

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        # change similarity matrix entries (i,j) to max{matrix}-(i,j)
        sim_max = numpy.max(numpy.ravel(sim_matrix))
        dist_matrix = sim_max - sim_matrix

        # set diagonal elements of the similarity matrix to zero
        for i in range(len(dist_matrix)):
            dist_matrix[i, i] = 0

        return dist_matrix

    def _guide_tree(self, dist_matrix):

        """
        @summary: Build a guide tree from the distance matrix

        @param dist_matrix: The distance matrix
        @type dist_matrix: numpy.ndarray
        @return: Pycluster similarity tree
        @rtype: Pycluster.cluster.Tree

        @author: Woon Wai Keen
        @author: Vladimir Likic
        """

        n = len(dist_matrix)

        print " -> Clustering %d pairwise alignments." % (n * (n-1)),
        tree = Pycluster.treecluster(distancematrix=dist_matrix, method='a')
        print "Done"

        return tree