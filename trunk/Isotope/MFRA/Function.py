"""
Functions for MFRA
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

from numpy import *

def c_mass_isotope_distr(mdv, c_corr):

    # calculate the exclusive mass isotope distribution of
    # the carbon skeleton, mdv_alpha_star
    mdv = array(mdv, float)
    mdv_alpha = mdv/sum(mdv)
    mdv_alpha = matrix(mdv_alpha)
    mdv_alpha = mdv_alpha.T
    c_corr = matrix(c_corr)
    mdv_alpha_star = (c_corr.I) * mdv_alpha
    mdv_alpha_star = array(mdv_alpha_star, float)
    mdv_alpha_star = mdv_alpha_star/sum(mdv_alpha_star)

    return mdv_alpha_star

def correction_matrix(n, num_a, nsil):

    """
    @summary: Calculates a correction matrix.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param num_a: Number of C, O, N, H, Si, or S atoms in the fragment (note:
        the number of C atoms excludes C atoms which may contain exogenous
        (non natural abundance) 13C).
    @type num_a: types.IntType
    @param nsil: Contains abundance values of natural stable isotopes.
    @type nsil: types.ListType
    
    @return: Correction matrix.
    @rtype: numpy.ndarray
    
    @author: Milica Ng
    """

    # if not(isinstance(n, IntType)):
        # error("'n' must be an integer")
    # elif not (n > 0):
        # error("'n' must be grater than zero")
    # if not(isinstance(num_a, IntType)):
        # error("'num_a' must be an integer")
    # elif not (num_a > 0):
        # error("'num_a' must be grater than zero")
    # if not(isinstance(nsil, ListType)):
        # error("'nsil' must be a list")
    # elif (nsil == []):
        # error("'nsil' must not be empty")
    # else:
        # for q in nsil: 
           # if not(isinstance(q, FloatType)):
               # error("'nsil' must be a list of decimal numbers")

    if not num_a == 0:
        c_corr_a = zeros((n+1,n+1))
        # calculate the first column of the correction matrix
        tmp = mass_dist_vector(n, num_a, nsil)
        c_corr_a[:,0]  = reshape(tmp,(1,len(tmp)))
        # calculate the rest of the matrix from the first column
        for i in range(1,(n+1)):
            for j in range(1,(n+1)):
                c_corr_a[i,j] = c_corr_a[i-1,j-1]
    else:
        c_corr_a = identity(n+1)
    return c_corr_a

def mass_dist_vector (n, num_a, nsil):

    """
    @summary: Calculates mass distribution vector.

    @param n: Number of fragment's C atoms which may contain exogenous (non
        natural abundance) 13C (e.g. only amino acid ones)
    @type n: types.IntType
    @param num_a: Number of C, O, N, H, Si, or S atoms in the fragment
        (note: the number of C atoms excludes C atoms which may contain
        exogenous (non natural abundance) 13C).
    @type num_a: types.IntType
    @param nsil: Contains abundance values of natural stable isotopes.
    @type nsil: types.ListType of types.FloatType

    @return: Mass distribution vector.
    @rtype: numpy.ndarray

    @author: Milica Ng
    """

    # if not(isinstance(n, IntType)):
        # error("'n' must be an integer")
    # elif not (n > 0):
        # error("'n' must be grater than zero")
    # if not(isinstance(num_a, IntType)):
        # error("'num_a' must be an integer")
    # elif not (num_a > 0):
        # error("'num_a' must be grater than zero")
    # if not(isinstance(nsil, ListType)):
        # error("'nsil' must be a list")
    # elif (nsil == []):
        # error("'nsil' must not be empty")
    # else:
        # for q in nsil: 
           # if not(isinstance(q, FloatType)):
               # error("'nsil' must be a list of decimal numbers")

    mdv = zeros(((n + 1),1))
    #   calculate all mass combinations where value 0 corresponds to m0,
    #     1 corresponds to m1, 2 corresponds to m2, ..., N corresponds to mN
    m_combs = []
    m_combs = combinations_with_repetition(size(range(0,(len(nsil)))), num_a) - 1
    # sort m_comb in order of row sums
    m_combs = m_combs.tolist()
    m_combs.sort(lambda x, y: cmp(sum(x),sum(y)))
    m_combs = array(m_combs)
    m = 0  # initial mass value
    count = 0  # current row in m_comb
    for i in range(0,n+1):
        #loop until the value of m changes (e.g m0 -> m1)
        while m == sum(m_combs[count,]):
            # calculate v vector
            v = zeros((len(nsil),1))
            for j in range(0,num_a):
                v[(m_combs[count,j])] = v[(m_combs[count,j])] + 1
            # calculate isotopolog abundance (part_1 * part_2) and add any
            # previous value if more than one combination corresponds to mN
            part_1 = fact(int(sum(v)))
            part_2 = 1
            for k in range(0,(len(nsil))):
                part_2 = part_2 * nsil[k]**v[k]/fact(int(v[k]))
            mdv[i] = mdv[i] + part_1 * part_2
            count = count + 1
            if (count > (size(m_combs,0)) - 1):  # check if finished
                break                              # (some entries might 0)
        if (((i + 1) > (size(m_combs,0)) - 1) or (count > (size(m_combs,0)) - 1)):
            break                              # (some entries might be 0)
        m = m + 1
    return mdv

def fact(integer):
    # from function import fact
    factorial = 1
    for factor in range(2, integer+1):
        factorial *= factor
    return factorial

def combinations_with_repetition(n, k):

    """
    @summary: Calculates indices for picking k elements (with replacement,
        order matters) from the set with n elements.

    @param n: Number of elements in the set
    @type n: types.IntType
    @param k: Number of picks
    @type k: types.IntType

    @return: Combinations with repetition indices (base 1).
    @rtype: numpy.ndarray

    @author: Milica Ng
    """

    # if not(isinstance(n, IntType)):
        # error("'n' must be an integer")
    # elif not (n > 0):
        # error("'n' must be grater than zero")
    # if not(isinstance(k, IntType)):
        # error("'k' must be an integer")
    # elif not (k > 0):
        # error("'k' must be grater than zero")

    if (k==1):
        indices = arange(1,n+1)[:, newaxis]
        return indices
    if (n==1):
        indices = ones((1,k,))
        return indices
    indices = []
    for z in range(1,n+1):
        next_indices = combinations_with_repetition(n+1-z,k-1)
        if indices == []:
            indices = hstack((z*ones((size(next_indices,0),1)), next_indices+z-1))
        else:
            tmp = hstack((z*ones((size(next_indices,0),1)), next_indices+z-1))
            indices = vstack((indices, tmp))
    return indices

