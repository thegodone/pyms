"""Array.py
 Module Array in metab.Utils
 Provides some convenience functions for lists/numpy
 arrays.
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

import math, types, copy

from Error import error
from Utils import *

def init_vec(N):

    """init_vec(N)
    Creates a vector of floats of dimension n initialized
    to 0.0, as a list.

    @param N An integer, vector dimension.
    @return Returns a list.
    """

    v = []
    for ii in range(N):
        v.append(0.0)
    return v

def init_mat(N, M):

    """init_mat(N, M)

    Creates an N x M matrix of floats initialized to 0.0,
    as a list of lists.

    @param N An integer. Number of rows
    @param M An integer. Dimension of each row 
    @return Returns a list.
    """

    m = []
    for ii in range(N):
        m.append(init_vec(M))
    return m

def transpose(m):

    """transpose(m)

    Transpose list of lists.
 
    @param m A list
    @return Returns a list of lists
    """

    if not is_list(m):
        error("argument not a list")
    else:
        for item in m:
            if not is_list(item):
                error("row not a list")

    M = len(m)
    N = len(m[0])
    mT = []

    for ii in range(N):
        rowT = []
        for jj in range(M):
            rowT.append(m[jj][ii])
        mT.append(rowT)

    return mT

def amean(v):

    """amean(v)

    Calculates mean.

    @param v A list or array
    @return A float
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    s = 0.0
    for e in v:
        s = s + e 
    s_mean = s/float(len(v))

    return s_mean

def astd(v):

    """astd(v)

    Calculates standard deviation.

    @param v A list or array
    @return A float
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    v_mean = amean(v)

    s = 0.0 
    for e in v:
        d = e - v_mean
        s = s + d*d
    s_mean = s/float(len(v)-1)
    v_std = math.sqrt(s_mean)

    return v_std

def amad(v):

    """amad(v)

    Calculates median absolute deviation.

    @parm v A list or array
    @return A float
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    m = amedian(v)
    m_list = []

    for xi in v:
        d = abs(xi - m)
        m_list.append(d)

    mad = amedian(m_list)/0.6745

    return mad

def amedian(v):

    """amediana(v)

    Calculates median.

    @parm v A list or array
    @return A float 
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    local_data = copy.deepcopy(v)
    local_data.sort()
    N = len(local_data)

    if (N % 2) == 0:
        # even number of points
        K = N/2 - 1 
        median = (local_data[K] + local_data[K+1])/2.0
    else:
	    # odd number of points
        K = (N - 1)/2 - 1
        median = local_data[K+1]

    return median

def amax(v):

    """amax(v)

    Finds the maximum element in a list or array.

    @param v A list or array
    @return Returns tuple (maxi, maxv), where maxv is the maximum
    element in the list and maxi is its index.
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    maxv = min(v) # Built-in min() function.
    maxi = None

    for ii in range(len(v)):
        if v[ii] > maxv:
            maxv = v[ii]
            maxi = ii

    if maxi == None:
        error("finding maximum failed")

    return maxi, maxv

def amin(v):

    """amin(v)

    Finds the minimum element in a list or array.

    @param v A list or array
    @return Returns tuple (maxi, maxv), where maxv is the minimum 
    element in the list and maxi is its index.
    """

    if not list_or_array(v):
        error("argument neither list nor array")

    minv = max(v) # Built-in max() function.
    mini = None

    for ii in range(len(v)):
        if v[ii] < minv:
            minv = v[ii]
            mini = ii

    if mini == None:
        error("finding maximum failed")

    return mini, minv

def costheta(v1, v2):

    """costheta(v1, v2)

    Returns cos theta for two vectors.

    @param v1 A list of numbers (vector 1)
    @param v2 A list of numbers (vector 2)
    @return A float, cosine theta between v1 and v2
    """

    if not is_list(v1):
        error("first argument is not a list")

    if not is_list(v2):
        error("second argument is not a list")

    if len(v1) != len(v2):
        error("vectors not of the same length")

    dp0 = dotp(v1, v2)
    dp1 = dotp(v1, v1)
    dp2 = dotp(v2, v2)
 
    ct = dp0/(math.sqrt(dp1)*math.sqrt(dp2))

    return ct 

def dotp(v1, v2):

    """dotp(v1, v2)

    Returns inner product of two vectors.

    @param v1 A list of numbers (vector 1)
    @param v2 A list of numbers (vector 2)
    @return A float, inner product of v1 and v2
    """

    if not is_list(v1):
        error("first argument is not a list")

    if not is_list(v2):
        error("second argument is not a list")

    if len(v1) != len(v2):
        error("vectors not of the same length")

    sumprod = 0.0

    for ii in range(len(v1)):
        sumprod = sumprod + v1[ii]*v2[ii]

    return sumprod

