"""
Time conversion functions
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

from Error import error
from Utils import is_str, is_str_num

def time_str_secs(time_str):

    """time_str_secs(time_str)

    Resolves time string of the form "<NUMBER>s" or "<NUMBER>m",
    returns time in seconds.

    time_str A string. Must be of the form "<NUMBER>s" or
        "<NUMBER>m" where "<NUMBER>" is a valid number.
    returns time in seconds
    """

    if not is_str(time_str):
        error("time string not a string")

    time_number = time_str[:-1]
    time_spec = time_str[-1].lower()

    if not is_str_num(time_number):
       print " --> received time string '%s'" % (time_number)
       error("improper time string")

    if not time_spec == "s" and not time_spec == "m":
        error("time string must end with either 's' or 'm'")

    time = float(time_number)

    if time_spec == "m":
        time = time*60.0

    return time

