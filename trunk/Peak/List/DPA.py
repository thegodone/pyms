"""DPA.py
"""

import numpy
import Pycluster

from pyms.Utils.Error import *

from pyms import Experiment

############################# Alignment class #################################

class Alignment:

    """ Holds a list of Alignments """

    def __init__(self, experiment=None):

        """ Initialization

        @param experiment Provide this if we're making a singleton alignment
            from an Experiment (optional)
        """

        if experiment and isinstance(experiment, Experiment.Class.Experiment):
            self.alignments = [experiment.peaks[:]] # work on copy of the peaklist
            self.expr_codes = [experiment.expr_code]
            self.similarity = 0.

    def __len__(self):
        return len(self.alignments)

    def write_csv(self, rt_filename, area_filename, minutes=True):

        """ Writes the alignments into a CSV file

        @param rt_filename File to save retention times to
        @param area_filename File to save normalized areas to
        @param minutes True/False to save retention time in minutes
        """

        try:
            fp_rt = open(rt_filename, 'w')
            fp_area = open(area_filename, 'w')
        except IOError:
            error("Cannot open output file for writing")

        # write experiment headers
        header = '"' + '","'.join(self.expr_codes) + "\"\n"

        fp_rt.write(header)
        fp_area.write(header)

        # for each peak, write each alignment's peak and area value
        for peak_idx in range(len(self.alignments[0])):
            rts = []
            areas = []
            for align_idx in range(len(self.alignments)):
                peak = self.alignments[align_idx][peak_idx]
                if peak is not None:
                    if minutes:
                        rtmin = peak.rt/60.0
                    else: 
                        rtmin = peak.rt
                    rts.append('%.3f' % rtmin)
                    areas.append('%.4f' % peak.norm_area)
                else:
                    rts.append('NA')
                    areas.append('NA')

            fp_rt.write(",".join(rts) + "\n")
            fp_area.write(",".join(areas) + "\n")

        fp_rt.close()
        fp_area.close()

    def transpose(self):

        """ Transposes the alignment lists (i.e. makes one row per peak match VS
            one row per experiment. For convenient processing by some routines. 
        """

        self.alignments = [ [r[col] for r in self.alignments]
                            for col in range(len(self.alignments[0]))
                          ]

    def filter_min_peaks(self, min_peaks):

        """ Filters away peak matches that have less than min_peaks

        @param min_peaks Minimum # of peaks
        """

        self.transpose()
        self.alignments =  [ x for x in self.alignments
                             if len(filter(None, x)) >= min_peaks
                           ]
        self.transpose()

################################ Tree class ###################################

class Tree:

    """ Holds a guide tree to be used with align_with_tree() """

    def __init__(self, reps, D, gap):
        
        """ Initialization

        @param reps A list of replicates. Can consist of Alignments or Experiements.
        @param D D parameter
        @param gap Gap parameter
        """

        self.D = D
        self.gap = gap

        # Handle trivial cases of 0 and 1
        if len(reps) == 0:
            return

        if len(reps) == 1:
            return

        n = len(reps)
        total = n * (n - 1) / 2

        # Step 1 - do pairwise alignments and make similarity matrix
        print " Calculating similarity matrix for %d replicates (D=%.2f, gap=%.2f)" % (n, D, gap)
        sim_matrix = numpy.zeros((n,n), dtype='f')
        for i in range(n - 1):
            for j in range(i + 1, n):
                alignment = align(reps[i], reps[j], D, gap)
                sim_matrix[i,j] = sim_matrix[j,i] = alignment.similarity

                total = total - 1
                print " -> %d replicate pairs remaining" % total

	# Step 2 - change similarity matrix entries (i,j) to: max {matrix} - (i,j)
	sim_max = numpy.max(numpy.ravel(sim_matrix))
	sim_matrix = sim_max - sim_matrix
	for i in range(n):
		sim_matrix[i,i] = 0

        self.sim_matrix = sim_matrix

        # Step 3 - perform hierchical clustering
        self.tree = Pycluster.treecluster(distancematrix=sim_matrix, method='a')

############################### public functions ##############################

def align_with_tree(As, T, min_peaks=1):

    """ Aligns a list of alignments/experiments As using a guide tree T.

    @return Merged alignments
    """

    print " Aligning %d replicates with guide tree (D=%.2f, gap=%.2f)" % (len(As), T.D, T.gap)

    # handle trivial cases of length 0 and 1
    if len(As) == 0:
        return []

    if len(As) == 1:
        return As

    # For everything else, we align according to the guide tree provided by
    # Pycluster. Quote documentation:
    #   Each item and subnode is represented by an integer. For hierarchical
    #   clustering of n items, we number the original items {0, ... , n-1},
    #   nodes are numbered {-1, ... , -(n-1)}. Note that the number of nodes
    #   is one less than the number of items.

    # extend As to length 2n to hold the n items, n-1 nodes, and 1 root
    As = As + [ None for _ in range(len(As)) ]

    # now align the alignments into positions -1, ... ,-(n-1)
    total = len(T.tree)
    index = 0
    for node in T.tree:
        index = index - 1
        As[index] = align(As[node.left], As[node.right], T.D, T.gap)

        total = total - 1
        print " -> %d replicates remaining" % total

    # our final alignment is in the root. filter min peaks and return..
    As[index].filter_min_peaks(min_peaks)
    return As[index]

def align(A1, A2, D, gap, sim_method="new"):

    """ Aligns two alignments/experiments

    @param A1 The first alignment/experiment
    @param A2 The second alignment/experiment
    @param D "D" parameter
    @param gap Gap parameters
    @param sim_method Use either "new" or "old" methods (default: new)

    @return Merged alignments
    """

    if (isinstance(A1, Experiment.Class.Experiment)):
        A1 = Alignment(A1)

    if (isinstance(A2, Experiment.Class.Experiment)):
        A2 = Alignment(A2)

    # calculate score matrix for these two runs
    M = score_matrix(A1, A2, D)

    # run dp algorithm on this matrix
    result = dp(M, gap)

    # make composite alignment from the results of the dp algo
    merged = merge_alignments(A1, A2, result['trace'])

    # calculate similarity score
    if sim_method == "old":
        cost = 0.
        for match in result['matches']:
            cost = cost + M[match[0]][match[1]]
        n = len(result['matches'])
        cost = cost/n * (len(A1.alignments[0]) + len(A2.alignments[0])) / (2*n)
        merged.similarity = cost
    else:
        merged.similarity = alignment_similarity(result['trace'], M, gap)

    return merged

########################### private functions #################################

def score_matrix(A1, A2, D):

    """ Calculates score matrix between two alignments. Called by align()

    @param A1 The first alignment/experiment
    @param A2 The second alignment/experiment
    @param D "D" parameter

    @return Merged alignments
    """

    score = numpy.zeros((len(A1.alignments[0]), len(A2.alignments[0])))
    times = numpy.zeros((len(A1.alignments[0]), len(A2.alignments[0])))

    # for all pairs of peaklists, calculate distance between them
    for peaks1 in A1.alignments:
        for peaks2 in A2.alignments:
            (score_matrix, binary_matrix) = score_matrix_pair(peaks1, peaks2, D)
            score = score + score_matrix
            times = times + binary_matrix

    # divide score_matrix by # of times the pair's score was calculated
    score = score / times

    return score

def score_matrix_pair(peaks1, peaks2, D):

    """ Calculates score matrix between two peaklists. Called by score_matrix()

    @param peaks1 The first list of peaks
    @param peaks2 The second list of peaks
    @param D "D" parameter

    @return A (score matrix, binary matrix) tuple
    """

    # Initialize gap matrix, Gap1 and Gap2. This is the identity matrix with
    # zero rows/columns inserted where there's a gap in peaks1/peaks2.
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
    binary_matrix = numpy.dot(numpy.dot(gap1, numpy.ones([len(f_peaks1), len(f_peaks2)])), gap2)

    return (score_matrix, binary_matrix)

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

    cs = 1. - (top / bot)

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

def dp(S, gap_penalty, peak_list1=None, peak_list2=None):
    
    """ Solves optimal path in score matrix based on global sequence alignment

    @param S Score matrix
    @param gap_penalty Gap penalty
    @param peak_list1 Peak list that coresponds to row index in score matrix
    @param peak_list2 Peak list that coresponds to column index in score matrix

    @return A dictionary of results
    """
    
    if(peak_list1 and peak_list2==None) or (peak_list1==None and peak_list2):
        error("Both peak lists must be defined")
    
    row_length=len(S[:,0])
    col_length=len(S[0,:])

    #D contains the score of the optimal alignment
    D = numpy.zeros((row_length+1,col_length+1), dtype='d')
    for i in range(1,row_length+1):
        D[i,0] = gap_penalty*i
    for j in range(1,col_length+1):
        D[0,j] = gap_penalty*j
    D[0,0] = 0.
    D[1:(row_length+1), 1:(col_length+1)] = S.copy();

    # Directions for trace
    # 0 - match               (move diagonal)
    # 1 - peaks1 has no match (move up)
    # 2 - peaks2 has no match (move left)
    # 3 - stop
    trace_matrix = numpy.zeros((row_length+1,col_length+1))
    trace_matrix[:,0] = 1; 
    trace_matrix[0,:] = 2;
    trace_matrix[0,0] = 3;
   
    for i in range(1,row_length+1):
        for j in range(1,col_length+1):
 
            #
            # Needleman-Wunsch Algorithm assuming a score function S(x,x)=0
            #
            #              | D[i-1,j-1] + S(i,j)
            # D[i,j] = min | D(i-1,j] + gap
            #              | D[i,j-1] + gap
            #

            darray = [D[i-1,j-1]+S[i-1,j-1], D[i-1,j]+gap_penalty, D[i,j-1]+gap_penalty]
            D[i,j] = min(darray)
            direction = darray.index(D[i,j])

            if(peak_list1 and peak_list2 and direction != 0):
                #Make sure gaps are inserted in order of retention time
                if(peak_list1[i-1].rt < peak_list2[j-1].rt):
                    direction = 2
                    D[i,j] = D[i,j-1]+gap_penalty
                elif(peak_list1[i-1].rt > peak_list2[j-1].rt):
                    direction = 1
                    D[i,j] = (D[i-1,j]+gap_penalty)

            trace_matrix[i,j] = direction

    # Trace back from bottom right
    trace = []
    matches = []
    i = row_length
    j = col_length
    direction = trace_matrix[i,j]
    p = [row_length-1]
    q = [col_length-1]
    
    while direction != 3:
        
        if direction == 0: #Match
            i=i-1
            j=j-1
            matches.append([i,j])
        elif direction == 1: #peaks1 has no match
            i=i-1
        elif direction == 2: #peaks2 has no match
            j=j-1
        p.append(i-1)
        q.append(j-1)
        trace.append(direction)
        direction=trace_matrix[i,j]

    #remove 'stop' entry
    p.pop()
    q.pop()
    # reverse the trace back
    p.reverse()
    q.reverse()
    trace.reverse()
    matches.reverse()
    return {'p':p, 'q':q, 'trace':trace, 'matches':matches, 'D':D, 'phi':trace_matrix}

def merge_alignments(A1, A2, traces):

    """ Merges two alignments with gaps added in from DP traceback

    @param A1 First alignment
    @param A2 Second alignment
    @param traces DP traceback

    @return A single alignment from A1 and A2
    """

    # Create object to hold new merged alignment and fill in its expr_codes
    merged_alignments = Alignment()
    merged_alignments.expr_codes = A1.expr_codes + A2.expr_codes

    # create empty lists of dimension |A1| + |A2|
    merged = [ [] for _ in range(len(A1) + len(A2)) ]

    A1 = A1.alignments
    A2 = A2.alignments

    idx1 = idx2 = 0
    for trace in traces:
        # trace can either be 0, 1, or 2
        # if it is 0, there are no gaps. otherwise, if it is 1 or 2,
        # there is a gap in A2 or A1 respectively.
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

    merged_alignments.alignments = merged

    # helper function for sorting alignments by its average
    def alignment_compare(x, y):
        x = [ _.rt for _ in filter(None, x)]
        y = [ _.rt for _ in filter(None, y)]

        avg_x = numpy.sum(x) / len(x)
        avg_y = numpy.sum(y) / len(y)

        if avg_x < avg_y:
            return -1
        else:
            return 1

    # sort according to average peak
    merged_alignments.transpose()
    merged_alignments.alignments.sort(alignment_compare)
    merged_alignments.transpose()

    return merged_alignments

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
    for trace in traces:
        # Trace can either be 0, 1, or 2
        # If it is 0, there is a match and we add to the sum the score between
	# these two aligned peaks. 
	# 
	# Otherwise, if it is 1 or 2, and there is a gap in A2 or A1
	# respectively. We then subtract the gap penalty from the sum.
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
