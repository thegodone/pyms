"""
Dynamic Programming routine
"""

import numpy

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

