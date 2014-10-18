"""
Project 4 - Computing alignments of sequences
Implement four functions:

build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score)
compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag)
compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)
compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)

"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Takes as input a set of characters alphabet and three scores diag_score, 
    off_diag_score, and dash_score. The function returns a dictionary of 
    dictionaries whose entries are indexed by pairs of characters in alphabet plus '-'.
    """

    alphabet_w_dash = set([])
    alphabet_w_dash = alphabet.copy()
    alphabet_w_dash.add('-')
   
    scoring_matrix = {}
    for row in alphabet_w_dash:
        scoring_matrix[row] = {}
        for idx in alphabet_w_dash:
            if row == idx:
                scoring_matrix[row][idx] = diag_score
            if row == '-' or  idx == '-':
                scoring_matrix[row][idx] = dash_score
            elif row != idx:
                scoring_matrix[row][idx] = off_diag_score

    return scoring_matrix

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common 
    alphabet with the scoring matrix scoring_matrix. The function computes and 
    returns the alignment matrix for seq_x and seq_y
    """

    alignment = [[0 for dummycol in range(len(seq_y)+1)] 
                           for dummyrow in range(len(seq_x)+1)]
    
    for idx in range(1, len(seq_x)+1):
        alignment[idx][0] = alignment[idx-1][0] + scoring_matrix[seq_x[idx-1]]['-']
        if global_flag == False and alignment[idx][0] < 0:
            alignment[idx][0] = 0
    for jdx in range(1, len(seq_y)+1):
        alignment[0][jdx] = alignment[0][jdx-1] + scoring_matrix['-'][seq_y[jdx-1]]
        if global_flag == False and alignment[0][jdx] < 0:
            alignment[0][jdx] = 0
    for idx in range(1, len(seq_x)+1):
        for jdx in range(1, len(seq_y)+1):
            alignment[idx][jdx] = max(
                                      alignment[idx-1][jdx-1] + scoring_matrix[seq_x[idx-1]][seq_y[jdx-1]],
                                      alignment[idx-1][jdx] + scoring_matrix[seq_x[idx-1]]['-'],
                                      alignment[idx][jdx-1] + scoring_matrix['-'][seq_y[jdx-1]])
            if global_flag == False and alignment[idx][jdx] < 0:
                alignment[idx][jdx] = 0
    
    return alignment

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet
    with the scoring matrix scoring_matrix. This function computes a global alignment 
    of seq_x and seq_y using the global alignment matrix alignment_matrix.
    """
    
    idx = len(seq_x)
    jdx = len(seq_y)
    seq_x_align = ''
    seq_y_align = ''
    
    while idx != 0 and jdx != 0:
        if alignment_matrix[idx][jdx] == alignment_matrix[idx-1][jdx-1] + scoring_matrix[seq_x[idx-1]][seq_y[jdx-1]]:
            seq_x_align = seq_x[idx-1] + seq_x_align
            seq_y_align = seq_y[jdx-1] + seq_y_align
            idx = idx - 1
            jdx = jdx - 1
        else:
            if alignment_matrix[idx][jdx] == alignment_matrix[idx-1][jdx] + scoring_matrix[seq_x[idx-1]]['-']:
                seq_x_align = seq_x[idx-1] + seq_x_align
                seq_y_align = '-' + seq_y_align
                idx = idx - 1
            else:
                seq_x_align = '-' + seq_x_align
                seq_y_align = seq_y[jdx-1] + seq_y_align
                jdx = jdx - 1
    while idx != 0:
        seq_x_align = seq_x[idx-1] + seq_x_align
        seq_y_align = '-' + seq_y_align
        idx = idx - 1
    while jdx != 0:
        seq_x_align = '-' + seq_x_align
        seq_y_align = seq_y[jdx-1] + seq_y_align
        jdx = jdx -1
    return (alignment_matrix[len(seq_x)][len(seq_y)], seq_x_align, seq_y_align)
               
def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common 
    alphabet with the scoring matrix scoring_matrix. This function computes a local 
    alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    """
    # find location of maximum value in alignment matrix
    max_value = 0
    x_coord, y_coord = (0,0)
    for idx in range(1, len(seq_x)+1):
        for jdx in range(1, len(seq_y)+1):
            if alignment_matrix[idx][jdx] > max_value:
                max_value = alignment_matrix[idx][jdx]
                x_coord, y_coord = idx, jdx
    
    idx = x_coord
    jdx = y_coord
    seq_x_align = ''
    seq_y_align = ''
    while alignment_matrix[idx][jdx] != 0:
        if alignment_matrix[idx][jdx] == alignment_matrix[idx-1][jdx-1] + scoring_matrix[seq_x[idx-1]][seq_y[jdx-1]]:
            seq_x_align = seq_x[idx-1] + seq_x_align
            seq_y_align = seq_y[jdx-1] + seq_y_align
            idx = idx - 1
            jdx = jdx - 1
        else:
            if alignment_matrix[idx][jdx] == alignment_matrix[idx-1][jdx] + scoring_matrix[seq_x[idx-1]]['-']:
                seq_x_align = seq_x[idx-1] + seq_x_align
                seq_y_align = '-' + seq_y_align
                idx = idx - 1
            else:
                seq_x_align = '-' + seq_x_align
                seq_y_align = seq_y[jdx-1] + seq_y_align
                jdx = jdx - 1
    
    return (max_value, seq_x_align, seq_y_align)         
    
def test():
    """
    Test functions
    """
    #testset = set(['a' ,'c', 't', 'g'])
    #scoring_matrix = build_scoring_matrix(testset, 10, 4, -4)
    #compute_alignment_matrix("act", "gta", scoring_matrix, True)
    #compute_global_alignment('A', 'A', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2}, '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, [[0, -4], [-4, 6]])
    #compute_local_alignment('', '', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2}, '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, [[0]])
test()