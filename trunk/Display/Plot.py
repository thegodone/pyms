"""Plot.py
 Module Plot in metab.Display.
 Provides plotting functions.
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

from metab.Utils.Utils import *
from metab import IO, Peak

import pylab
import numpy

__hover_line = None  # Line2D instance of line that follows mouse over graph
__current_line = None  # Line2D instance of red line showing rt of spectral view
__figure_data = None  # ChemStation instance for current graph

def ic(ics, labels=None, minutes=False, markers=False):

    """ic(ics)

    Plot a list of IonChromatogram instances.

    @param ics A single or a list of IonChromatograms for plotting
    @param labels A list of line labels (default: None)
    @param minutes Set to True to show x-axis in minutes (default: show
        in seconds)
    @param ticks Set to True to show markers at scan points
    """

    # If we don't get a list to the input to ics, we either handle it if it's
    # a IonChromatogram instance, or reject if it isn't.
    if not isinstance(ics, list):
        if isinstance(ics, IO.Class.IonChromatogram):
            ics = [ics]
        else:
            error("ics argument must be an IonChromatogram or a list of it")
        
    # Initialize graph
    pylab.grid(True)
    pylab.hold(True)
    if not labels == None:
        pylab.legend()
    __xaxisminutes(minutes)

    # Plot ICs
    for i in range(len(ics)):
        if not isinstance(ics[i], IO.Class.IonChromatogram):
            error("ics argument must be an IonChromatogram or a list of it")

        x = ics[i].get_time_array()
        y = ics[i].get_intensity_array()

        if not labels == None:
            pylab.plot(x, y, label=labels[i])
        else:
            pylab.plot(x, y)

        if markers:
            pylab.plot(x, y, 'rx')

    pylab.show()

def sm(data, minutes=False):

    """tic(data, minutes=False)

    Plot an instance of ChemStation and allows for inspection of mass spectrum

    @param data A ChemStation data structure
    @param minutes Set to True to show x-axis in minutes (default: show
        in seconds)
    """

    global __hover_line, __current_line, __figure_data

    # Type checking
    if not isinstance(data, IO.ANDI.ChemStation):
        error("data argument must be an instance of ChemStation")

    # Initialize graph
    pylab.figure(1)
    pylab.hold(True)
    pylab.grid(True)
    __xaxisminutes(minutes)

    # Get TIC from data and plot in black
    ic = data.get_tic();
    x = ic.get_time_array()
    y = ic.get_intensity_array()
    pylab.plot(x, y, 'k')

    # Set up handler for viewing m/z
    pylab.connect('button_release_event', __mouseclick)
    pylab.connect('motion_notify_event', __mousemove)
    __figure_data = data

    pylab.show()

    # Clean up globals
    __hover_line = None
    __current_line = None
    __figure_data = None

def peaks(ic, peaks, minutes=False):

    """peaks(ic, peaks, minutes=False)

    Plot an instance of IonChromatogram and its peaks.

    @param ic An IonChromatogram for plotting
    @param peaks A list of Peaks corresponding to the ion chromatogram
    @param minutes Set to True to show x-axis in minutes (default: show
        in seconds)
    """

    global __hover_line, __current_line, __figure_data

    # Input checking
    if not isinstance(ic, IO.Class.IonChromatogram):
        error("ic argument must be an IonChromatogram")

    # Initialize graph
    pylab.figure(1)
    pylab.hold(True)
    pylab.grid(True)
    __xaxisminutes(minutes)

    # Initialize local structures
    x = ic.get_time_array()
    y = ic.get_intensity_array()
    fillcolor = None
    p_xs = []
    p_ys = []

    # Prepare peak plots and peak area fill
    for peak in peaks:
        if not isinstance(peak, Peak.Class.Peak):
            error("peak argument must only consist of Peak objects")

        p_xs.append(peak.rt)
        p_ys.append(peak.intensity)
        p_left = peak.pt_bounds[0]
        p_peak = peak.pt_bounds[1]
        p_right = peak.pt_bounds[2]
    
        # Initialize polygon to lower-left corner, populate polygon array,
        # close the polygon, then fill it with alternate colors.
        p_x = [x[p_left]]
        p_y = [0]
    
        for i in range(p_left, p_right+1):
            p_x.append(x[i])
            p_y.append(y[i])

        p_x.append(x[p_right])
        p_y.append(0)

        if fillcolor == '0.8':
            fillcolor = '0.5'
        else:
            fillcolor = '0.8'
      
        pylab.fill(p_x, p_y, fillcolor)

        if minutes:
            p_label = '%s' % time_secs_to_mins(peak.rt)
        else:
            p_label = '%.2f' % peak.rt
        pylab.text(peak.rt, peak.intensity, p_label,
            horizontalalignment='center')
    
    # Plot peaks in red crosses and IC in black
    pylab.plot(p_xs, p_ys, 'rx')
    pylab.plot(x, y, 'k')

    pylab.show()

def mz(data, time, top=5):

    """mz(data, time)

    Plot mass spectrum at a given time from Chemstation data.

    @param data A ChemStation data structure
    @param time Time
    @param top Draw numbers for "top" number of peaks (default: 5)
    """

    # Type checking
    if not isinstance(data, IO.ANDI.ChemStation):
        error("data argument must be an instance of ChemStation")
    if not is_number(time):
        error("time argument must be numeric")

    # Initialize graph
    pylab.figure(2)
    pylab.hold(False)

    idx = data.get_index_at_time(time)
    mzs = data.get_mass_spectrum_at_index(idx)
    mz_idx_params = data.get_mass_range()
    mz_idxs = range(mz_idx_params[0], mz_idx_params[1]+1)
    
    pylab.bar(mz_idxs, mzs, 0.2)
    pylab.xlabel('m/z at time %.1f' % time)
    
    # Draw top n peaks by reverse sorting them by intensity
    sorted_mzs = zip(mzs, mz_idxs)
    sorted_mzs.sort()
    for mz in sorted_mzs[-top:]:
        (mass, idx) = mz
        pylab.text(idx, mass, str(idx), horizontalalignment='center')

    pylab.draw()
    pylab.show()

def __xaxisminutes(minutes=False):

    """__xaxisminutes(minutes=False)

    Private function for setting x-axis on current figure to minutes or seconds
    
    @param minutes Set to True to show x-axis in minutes (default: show
        in seconds)
    """

    if minutes:
        formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x, y:
                time_secs_to_mins(x))
    else:
        formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x, y:
                '%.1f' % x)
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(formatter)

def __mouseclick(event):

    """__mouseclick(event)

    Private callback function for drawing mass spectrum when TIC graph
    is clicked.
    
    @param event MplEvent object representing the event
    """

    global __figure_data, __current_line

    # On middle-click, draw mass spectrum at the point, and also a line on IC
    # indicating time of spectrum view
    if event.button == 2:
        mz(__figure_data, event.xdata)

        pylab.figure(1)
        if __current_line:
            __current_line.set_xdata(event.xdata)
        else:
            __current_line = pylab.axvline(event.xdata, color='r')

def __mousemove(event):

    """__mousemove(event)

    Private callback function for drawing line selection when mouse is
    moved over TIC graph.
    
    @param event MplEvent object representing the event
    """

    global __hover_line

    # Update line if mouse pointer is moving over the graph.
    if event.inaxes:
        pylab.figure(1)

        if __hover_line:
            __hover_line.set_xdata(event.xdata)
        else:
            __hover_line = pylab.axvline(event.xdata)

        pylab.draw()

