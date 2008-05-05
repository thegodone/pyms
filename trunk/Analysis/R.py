"""R.py
 Module R in metab.Utils
 Provides link to R statistical package.
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

from rpy import *

from metab.Utils.Error import error

def t_two_sample(a, b, equal_variance=True):

    """t_two_sample(a, b, equal_variance=True)

    Two-sample t-test to test the hypothesis that two samples come
    from the distribution from the same mean. Normal distrubution
    assumed.

    @param a A list of numbers.
    @param b A list of numbers.
    @return A list [ p-value, t, df ].
    """

    if equal_variance:
        R_logical = r.TRUE
    else:
        R_logical = r.FALSE

    t = r.t_test(a, b, var_equal=R_logical)

    q = [ t['p.value'], t['statistic']['t'], t['parameter']['df'] ]

    return q 

def wilcox_two_sample(a, b):

    """wilcox_two_sample(a, b)

    Two sample Wilcoxon rank sum test (also known as Mann-Whitney test).
    A non-parametric test used to compare two independent groups of
    sampled data. No assumptions about the distribution of the data
    is made.

    @param a A list of numbers.
    @param b A list of numbers.
    @return A list [ p-value, W-statistic ].
    """

    t = r.wilcox_test(a, b)
    q = [ t['p.value'], t['statistic']['W'] ]

    return q

