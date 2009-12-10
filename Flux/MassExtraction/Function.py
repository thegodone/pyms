"""
Functions for Flux.MassExtract
"""


 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-9 Vladimir Likic                                    #
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


from pyms.Flux.Class import MIDS
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

# for plotting only, delete later!
from pylab import * # used by plotfile, savefig

def peak_apex(ic, rt, win_size, compound, file_name, mz):

    """
    @summary: Locate peak apex, as a local maximum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param rt: Retention time
    @type rt: FloatType
    @param win_size: Window for a local maximum search
    @type win_size: FloatType
    @param compound: Compound name used to raise the warning 
    @type compound:StringType
    @param file_name: File name used to raise the warning
    @type file_name: TypeString
    @param mz: Mass to charge ratio used to raise the warning
    @type mz: IntType
    
    @return: Ion Chromatogram index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # get indices at retention time +- 1/2 window size
    lb_rt = rt - (win_size/2)
    lb_ix = ic.get_index_at_time(lb_rt)
    ub_rt = rt + (win_size/2)
    ub_ix = ic.get_index_at_time(ub_rt)

    # find peak_apex (local maximum between lb & ub indices)
    peak_apex_int = ic.get_intensity_at_index(lb_ix)
    peak_apex_ix = lb_ix

    for ix in range(lb_ix,ub_ix+1):
       if ic.get_intensity_at_index(ix) > peak_apex_int:
          peak_apex_int = ic.get_intensity_at_index(ix)
          peak_apex_ix = ix

    # peak_apex_rt = ic.get_time_at_index(peak_apex_ix)

    return peak_apex_ix

def left_boundary(ic, peak_apex_ix, peak_apex_int):

    """
    @summary: Locate left boundary, as a local minimum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param peak_apex_ix: Peak apex index
    @type peak_apex_ix: IntType
    @param peak_apex_int: Peak apex intensity
    @type peak_apex_ix: FloatType
    
    @return: Left boundary index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # find left_boundary (first local minimum to the left of peak apex)
    left_boundary_int = peak_apex_int
    ix = peak_apex_ix
    while ix >= 0:
        if ic.get_intensity_at_index(ix) > left_boundary_int:
            break
        else:
            left_boundary_int = ic.get_intensity_at_index(ix)
            ix = ix - 1
    left_boundary_ix = ix + 1
    # left_boundary_rt = ic.get_time_at_index(left_boundary_ix)   

    return left_boundary_ix

def right_boundary(ic, peak_apex_ix, peak_apex_int):

    """
    @summary: Locate right boundary, as a local minimum 

    @param ic: Ion Chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param peak_apex_ix: Peak apex index
    @type peak_apex_ix: IntType
    @param peak_apex_int: Peak apex intensity
    @type peak_apex_ix: FloatType
    
    @return: Right boundary index 
    @rtype: IntType
    
    @author: Milica Ng
    """

    # find right_boundary (first local minimum to the right of peak apex)
    right_boundary_int = peak_apex_int
    ix = peak_apex_ix
    while ix < len(ic):
        if ic.get_intensity_at_index(ix) > right_boundary_int:
            break
        else:
            right_boundary_int = ic.get_intensity_at_index(ix)
            ix = ix + 1
    right_boundary_ix = ix - 1
    # right_boundary_rt = ic.get_time_at_index(right_boundary_ix)
    return right_boundary_ix


def calculate_mid(ic_list, ave_left_ix, ave_right_ix):

    """
    @summary: Calculate mass isotopomer distribution

    @param ic_list: List of IonChromatograms
    @type ic_list: ListType
    @param ave_left_ix: Average left boundary index for the ions
    @type ave_left_ix: IntType
    @param ave_right_ix: Average right boundary index for the ions
    @type ave_right_ix: IntType
    
    @return: None
    @rtype: NoneType
    
    @author: Milica Ng
    """

    # calculate mass isotopomer distribution
    mid = []
    for ic in ic_list:
        # get ion chromatogram at specified m/z
        # ic = im.get_ic_at_mass(mz)
        area = 0
        # integrate (sum intensity from average left_boundary_ix to average right_boundary_ix 
        for ix in range(ave_left_ix, ave_right_ix+1):
            area = area + ic.get_intensity_at_index(ix)
        mid.append(area)

    return mid

def plot_ics(ic_list, ave_left_ix, ave_right_ix, file_name, compound, noise):

    """
    @summary: For testing purposes only! To be deleted.
    
    @author: Milica Ng
    """

    count = 0
    for ic in ic_list:
        col1 = []
        col2 = []

        for ix in range(ave_left_ix-11,ave_right_ix+11):
            col1.append(ic.get_time_at_index(ix))
            col2.append(ic.get_intensity_at_index(ix))

        plot_file = compound+"-"+str(file_name)+"-"+str(count)+".png"
        count = count+1        
        # print '-> Plotting ', plot_file, '\n'
        plot(col1,col2) 
        vlines(ic.get_time_at_index(ave_left_ix), 0, noise)
        vlines(ic.get_time_at_index(ave_right_ix), 0, noise)
        title(plot_file)
        savefig(plot_file)
        close()


def extract_mid(file_name, im, mids, win_size, noise):

    """
    @summary: Method for extracting mass isotopomer distribution (MID)

    @param file_name: File number
    @type file_name: StringType
    @param im: Intensity matrix
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param mids: Mass isotopomer distribution data
    @type mids: pyms.Flux.Class.MIDS
    @param mid_size: Size of the mass distribution vector
    @type mid_size: IntType
    @param win_size: Size of the window for local maximum search
    @type win_size: FloatType
    @param noise: Noise threshold below which signal is assumed unreliable
    @type win_size: IntType

    @return: Mass isotopomer distribution data
    @rtype: pyms.Flux.Class.MIDS
    
    @author: Milica Ng
    """
    
    # get name, rt, ion and MID size
    compound = mids.get_name()
    rt = mids.get_rt()
    ion = mids.get_ion()
    mid_size = mids.get_mid_size()

    # initialise loop variables
    left_boundary_sum = 0
    right_boundary_sum = 0
    sum_count = 0
    ic_list = []
    peak_apex_time_list = []
    peak_apex_intensity_list = []

    # loop over required number of mass isotopomers
    for mz in range(ion, ion+mid_size):

        # get ion chromatogram at current m/z
        # todo: check if mass is out of range, raise an error!
        ic = im.get_ic_at_mass(mz) 

        # smooth noise + correct baseline
        ic_smooth = savitzky_golay(ic)
        ic_bc = tophat(ic_smooth, struct="1.5m")
        ic = ic_bc

        # store ic's for calculate_mid later
        ic_list.append(ic)

        # get peak apex index and intensity
        peak_apex_ix = peak_apex(ic, rt, win_size,compound, file_name, mz)
        peak_apex_int = ic.get_intensity_at_index(peak_apex_ix)

        # store peak apex retention time & intensity for later checks
        peak_apex_time_list.append(ic.get_time_at_index(peak_apex_ix))
        peak_apex_intensity_list.append(ic.get_intensity_at_index(peak_apex_ix))

        # get indices for left/right boundary
        left_boundary_ix = left_boundary(ic, peak_apex_ix, peak_apex_int)
        right_boundary_ix = right_boundary(ic, peak_apex_ix, peak_apex_int)

        if peak_apex_int > noise: 

            # raise warning if current peak apex retention time is within 2secs of previous ones
            for time in peak_apex_time_list:
                if (abs(ic.get_time_at_index(peak_apex_ix) - time)) > 2:
                     print 'WARNING:', compound, ion,' The peak apex retention time shift is larger than 2 secs in' , file_name, 'for ion', mz
                break
            # raise warning if peak apex is smaller than its neighbours
            if ic.get_intensity_at_index(peak_apex_ix) < ic.get_intensity_at_index(peak_apex_ix-1):
                print 'WARNING:', compound, ion,' There is a larger intensity to the left of peak apex in' , file_name, 'ion', mz
            if ic.get_intensity_at_index(peak_apex_ix) < ic.get_intensity_at_index(peak_apex_ix+1):
                print 'WARNING:', compound, ion,'  There is a larger intensity to the right of peak apex in' , file_name, 'ion', mz

            # sum left and right boundary indices above noise in the current file
            left_boundary_sum = left_boundary_sum + left_boundary_ix
            right_boundary_sum = right_boundary_sum + right_boundary_ix
            sum_count = sum_count + 1

    if sum_count > 0:

        # raise warning if peak boundaries belonging to a single peak were used as average
        if sum_count == 1:
        print 'WARNING:',compound, ion,' Only one peak was larger than', noise, '!','in' , file_name

        # find average of left and right boundary index in the current file
        ave_left_ix = left_boundary_sum/sum_count
        ave_right_ix = right_boundary_sum/sum_count

        # raise warning if an intensity larger than noise is found where peak apex is smaller than noise
        j = 0
        for peak_apex_intensity in peak_apex_intensity_list:
             if peak_apex_intensity < noise:
                 for ix in range(ave_left_ix, ave_right_ix):
                     if ic_list[j].get_intensity_at_index(ix) > noise:
                         print 'WARNING:', compound,' An intensity greater than', noise, 'in' , file_name, 'ion',ion,'+', j
                         break
             j = j+1

        # plot for testing purposes only. to be deleted!
        plot_ics(ic_list, ave_left_ix, ave_right_ix, file_name, compound, noise)

        # calculate mass isotopomer distribution in the current file
        mid = calculate_mid(ic_list, ave_left_ix, ave_right_ix)

        # store mass isotopomer distribution inside mids variable
        mids.set_values(mid, file_name)
            
    else:

        mid = [0] * mid_size

        # store mass isotopomer distribution inside mids variable
        mids.set_values(mid, file_name)

    return mids
