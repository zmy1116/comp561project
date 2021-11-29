"""
This file contains an implementation of the Needleman-Wunsch algorithm with
affine gap penalty.
"""

# %% Useful imports

import numpy as np
from src.utils.utils import sequence_one_hot
from src.blast.nt_scoring_function import NT_SCORE_ALGORITHM


# %% Methods to decide what to put in one matrix slot


def get_max_M(M, X, Y, i, j):
    origins = {'X': X[(i, j)],
               'Y': Y[(i, j)],
               'M': M[(i, j)]}
    score = max(origins.values())
    res = [ind for ind in origins if origins[ind] == score]
    return (score, [res[0]])


def get_max_X(M, X, Y, i, j, gap_create, gap_extend, N):
    # Only extend if we are at the Nth nuc
    if (i+1 % N) == 0:
        origins = {'X': X[(i, j)] + gap_extend,
                   # No gaps in X if there was a gap in Y
                   #'Y': Y[(i, j)] + gap_create + gap_extend,
                   'M': M[(i, j)] + gap_create + gap_extend}
    else:
        origins = {'X': X[(i, j)] + gap_extend}
    score = max(origins.values())
    res = [ind for ind in origins if origins[ind] == score]
    return (score, [res[0]])


def get_max_Y(M, X, Y, i, j, gap_create, gap_extend, N):
    # Only extend if we are at the Nth nuc
    if (j+1 %N) == 0:
        origins = {# No gaps in Y if there was a gap in X
                   #'X': X[(i, j)] + gap_create + gap_extend,
                   'Y': Y[(i, j)] + gap_extend,
                   'M': M[(i, j)] + gap_create + gap_extend}
    else:
        origins = {'Y': Y[(i, j)] + gap_extend}
    score = max(origins.values())
    res = [ind for ind in origins if origins[ind] == score]
    return (score, [res[0]])


# %% Method to get the paths followed to compute a given matrix


def get_all_paths(lastConsidered, paths, draftPath, M_Previous, X_Previous,
                  Y_Previous):
    """
    Input:
    Output: all optimal paths leading to the bottom right square of the matrix
    """
    if lastConsidered[0] == (0, 0):
        paths.append(draftPath)
        return
    else:
        (i, j) = lastConsidered[0]
        if lastConsidered[1] == 'X':
            draftPathBis = draftPath + ['X']
            for origin in X_Previous[(i, j)]:
                lastCoords = (i, j - 1)
                lastConsideredBis = (lastCoords, origin)
                get_all_paths(lastConsideredBis, paths, draftPathBis,
                              M_Previous, X_Previous, Y_Previous)
        if lastConsidered[1] == 'Y':
            draftPathBis = draftPath + ['Y']
            for origin in Y_Previous[(i, j)]:
                lastCoords = (i - 1, j)
                lastConsideredBis = (lastCoords, origin)
                get_all_paths(lastConsideredBis, paths, draftPathBis,
                              M_Previous, X_Previous, Y_Previous)
        if lastConsidered[1] == 'M':
            draftPathBis = draftPath + ['M']
            if (i, j) not in M_Previous:
                return 'M'
            for origin in M_Previous[(i, j)]:
                lastCoords = (i - 1, j - 1)
                lastConsideredBis = (lastCoords, origin)
                get_all_paths(lastConsideredBis, paths, draftPathBis,
                              M_Previous, X_Previous, Y_Previous)


# %% Method to get one path followed to compute a given matrix

def get_one_path(lastConsidered, draftPath, M_Previous, X_Previous,
                 Y_Previous):
    """
    Input:
    Output: one optimal path leading to the bottom right square of the matrix
    """
    while lastConsidered[0] != (0, 0):
        (i, j) = lastConsidered[0]

        if lastConsidered[1] == 'X':
            draftPath = draftPath + ['X']
            origin = X_Previous[(i, j)][0]
            lastCoords = (i, j - 1)
            lastConsidered = (lastCoords, origin)

        if lastConsidered[1] == 'Y':
            draftPath = draftPath + ['Y']
            origin = Y_Previous[(i, j)][0]
            lastCoords = (i - 1, j)
            lastConsidered = (lastCoords, origin)

        if lastConsidered[1] == 'M':
            draftPath = draftPath + ['M']
            origin = M_Previous[(i, j)][0]
            lastCoords = (i - 1, j - 1)
            lastConsidered = (lastCoords, origin)

    return draftPath


# %% Matrix computation


def nw_affine_matrix_two(seq1, seq2, gap_create, gap_extend,
                         score_method, threshold_score, N,
                         mismatch_score=5,
                         substitution_dict=dict()):
    if threshold_score is None:
        threshold_score = np.infty
    """
    Input: sequences to align and alignment info
    Output: best possible alignments and corresponding score
    """
    p = len(seq1)
    n = len(seq2)
    M_Matrix = np.zeros((n + 1, p + 1))
    X_Matrix = np.zeros((n + 1, p + 1))
    Y_Matrix = np.zeros((n + 1, p + 1))
    M_Previous = dict()
    X_Previous = dict()
    Y_Previous = dict()

    seq1_onehot = sequence_one_hot(seq1)

    # Initialize first matrix cells (first column for X, first line for Y)
    for i in range(n + 1):
        X_Matrix[(i, 0)] = - np.infty
        if (i > 0):
            M_Matrix[(i, 0)] = - np.infty
            if (i == 1):
                Y_Matrix[(i, 0)] = gap_create + gap_extend
                Y_Previous[(i, 0)] = 'M'
            else:
                Y_Matrix[(i, 0)] = Y_Matrix[(i - 1, 0)] + gap_extend
                Y_Previous[(i, 0)] = 'Y'
    for j in range(p + 1):
        Y_Matrix[(0, j)] = - np.infty
        if (j > 0):
            M_Matrix[(0, j)] = - np.infty
            if (j == 1):
                X_Matrix[(0, j)] = gap_create + gap_extend
                X_Previous[(0, j)] = 'M'
            else:
                X_Matrix[(0, j)] = X_Matrix[(0, j - 1)] + gap_extend
                X_Previous[(0, j)] = 'X'

    # Now fill in the rest of the matrixes, recording the paths
    for i in range(n):
        letter2_probas = seq2[i]
        for j in range(p):
            letter1 = seq1_onehot[j]
            substitution_score = NT_SCORE_ALGORITHM[score_method](letter2_probas,
                                                                  letter1,
                                                                  mismatch_score,
                                                                  substitution_dict)
            # M
            M_Previous[(i + 1, j + 1)] = []
            score, resM = get_max_M(M_Matrix, X_Matrix, Y_Matrix, i, j)
            for origin in resM:
                M_Previous[(i + 1, j + 1)].append(origin)
                M_Matrix[(i + 1, j + 1)] = score + substitution_score
            # X
            X_Previous[(i + 1, j + 1)] = []
            score, resX = get_max_X(M_Matrix, X_Matrix, Y_Matrix, i + 1, j,
                                    gap_create, gap_extend, N)
            for origin in resX:
                X_Previous[(i + 1, j + 1)].append(origin)
                X_Matrix[(i + 1, j + 1)] = score
            # Y
            Y_Previous[(i + 1, j + 1)] = []
            score, resY = get_max_Y(M_Matrix, X_Matrix, Y_Matrix, i, j + 1,
                                    gap_create, gap_extend, N)
            for origin in resY:
                Y_Previous[(i + 1, j + 1)].append(origin)
                Y_Matrix[(i + 1, j + 1)] = score
            
            # Check if we need to stop
            curr_score = max([M_Matrix[(i + 1, j + 1)], X_Matrix[(i + 1, j + 1)], Y_Matrix[(i + 1, j + 1)]])
            if curr_score >= threshold_score:
                break
        if curr_score >= threshold_score:
            break
    
    
    final_score, final_mat = get_max_M(M_Matrix, X_Matrix, Y_Matrix, 0, p)
    final_pos = (0, p)

    for k in range(1, i + 1):
        tmp_score, tmp_mat = get_max_M(M_Matrix, X_Matrix, Y_Matrix, k, p)
        if tmp_score > final_score:
            final_score = tmp_score
            final_mat = tmp_mat
            final_pos = (k, p)

    return (M_Matrix, X_Matrix, Y_Matrix, M_Previous, X_Previous, Y_Previous,
            final_score, final_pos, final_mat)


# %% Method to get the sequence alignment corresponding to a given path


def from_path_to_als(path, seq1, seq2, offset):
    al1 = ""
    al2 = ""
    ind1 = len(seq1) - 1
    ind2 = len(seq2) - 1
    for position in path:
        if position == "M":
            al1 = seq1[ind1] + al1
            al2 = "ACGT"[np.argmax(seq2[ind2])] + al2
            ind1 -= 1
            ind2 -= 1
        elif position == "X":
            al1 = seq1[ind1] + al1
            ind1 -= 1
            al2 = "-" + al2
            if (ind2 >= 0):
                offset["seq2"] = offset["seq2"][0: ind2 + 1] + [-1] + offset["seq2"][ind2 + 1: len(offset["seq2"])]
            else:
                offset["seq2"] = [-1] + offset["seq2"]
        else:
            al2 = "ACGT"[np.argmax(seq2[ind2])] + al2
            ind2 -= 1
            al1 = "-" + al1
            if (ind1 >= 0):
                offset[seq1] = offset[seq1][0: ind1 + 1] + [-1] + offset[seq1][ind1 + 1: len(offset[seq1])]
            else:
                offset[seq1] = [-1] + offset[seq1]
    return (al1, al2)


# %%


def nw_affine_two(seq1, seq2, gap_create, gap_extend, score_method,
                  threshold_score, N,
                  substitution_dict=dict(), offset=0, all_paths=False,
                  mismatch_score=5):
    M_Matrix, X_Matrix, Y_Matrix, M_Previous, X_Previous, Y_Previous, finalScore, finalPos, finalMat = nw_affine_matrix_two(
        seq1,
        seq2,
        gap_create,
        gap_extend,
        score_method,
        threshold_score,
        N,
        mismatch_score,
        substitution_dict)
    paths = []
    p = len(seq1)
    n = len(seq2)
    if offset == 0:
        offset = dict()
        offset[seq1] = list(range(len(seq1)))
        offset["seq2"] = list(range(len(seq2)))
    for mat in finalMat:
        get_all_paths((finalPos, mat), paths, [], M_Previous, X_Previous,
                      Y_Previous)
    alignments = from_path_to_als(paths[0], seq1, seq2, offset)
    return (alignments, finalScore)


# %%

query = "AAT"
ref = "AAA"

gap_create = 0
gap_extend = -1

reference = np.zeros((3, 4))
reference[(0, 0)] = 1
reference[(1, 0)] = 1
reference[(2, 0)] = 1

mismatch_score = 10

# expect: AAT VS AA-A ou AA-T VS AAA

score_method = "sum_proba_score"
(alignments, finalScore) = nw_affine_two(query, reference, gap_create, gap_extend, score_method,
                                         N=1, threshold_score=-np.infty,
                                         substitution_dict=dict(), offset=0, all_paths=False,
                                         mismatch_score=5)

# %%
"""
reference = reference_matrix[14069:14147]
query = "A" * 26
s_query = query[12:]
score_method = "sum_proba_score"
nw_affine_two(s_query, reference, -1, -2, score_method,
              substitution_dict=dict(),offset=0, all_paths=False,
              mismatch_score=5)
"""

"""
pos_left_q = 15
query = "T" + "A" * 26
seq1 = query[:pos_left_q][::-1]
print(seq1)
max_length = 3 * len(query)
pos_left_ref = 593383
seq2 = reference_matrix[max(0, pos_left_ref - 1 - max_length):pos_left_ref][::-1]
print(seq2)

gap_create = -5
gap_extend = -1
mismatch_score = 5
score_method = "sum_proba_score"


(alignments, finalScore)=nw_affine_two(seq1, seq2, gap_create, gap_extend,
                                       score_method, substitution_dict=dict(),
                                       offset=0, all_paths=False,
                                       mismatch_score=5)
"""
