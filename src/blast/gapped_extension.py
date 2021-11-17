from src.blast import ungapped_extension
import pickle
import numpy as np
from src.utils.utils import sequence_one_hot

def gapped_alignment(query, reference, pos_q, pos_r, score_method, substitution,
                     mismatch_score, gap_penalty, gap_bias = 0, reverse = False):
    """
    compute a gapped alignment
    generates an alignment between the query and the reference sequences
    
    :param query: query sequence of letters (ACGT)
    :param reference: reference matrix file
    :param pos_q: position of the query at which to start the alignment
    :param pos_r: position of the reference at which to start the alignment
    :param score_method: method used to compute extension score
    :param substitution_dict: stores the cost of replacing one letter by another
    :param gap_penalty: cost of a gap in the alignment
    :param gap_bias: cost for opening a gap, if affine function
    :param reverse: if the alignment going downstream? (default: False)
    
    :return: alignment strings, the position on the reference, and the score
    """
    # reference_matrix = pickle.load(open(reference, 'rb'))
    reference_matrix = reference
    #Subset the query and reference matrix to get the parts to align
    #To reduce time, the length subsetted to the ref matrix is max 3 times the length of the sequence subset
    max_length = len(query)*3
    
        
    if reverse :
        s_query = query[:pos_q]
        start = np.max([0, pos_r-max_length])
        s_reference_matrix = reference_matrix[start:pos_r]
        # Read the string and matrix backwards
        s_query = s_query[::-1]
        s_reference_matrix = np.flipud(s_reference_matrix)
    else :
        s_query = query[pos_q:]
        end = np.min([pos_r+max_length, len(reference_matrix)])
        s_reference_matrix = reference_matrix[pos_r:end]
        
    s_query_one_hot = sequence_one_hot(s_query)
    sM = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    sX = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    sY = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pM = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pX = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pY = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    
    # Initialization
    sM[:,0] = - np.infty
    sM[0,:] = - np.infty
    sX[:,0] = - np.infty
    sY[0,:] = - np.infty
    sY[:,0] = [0] + [gap_bias + k * gap_penalty for k in range(1, len(s_query)+1)]
    sX[0,:] = [0] + [gap_bias + k * gap_penalty for k in range(1, len(s_reference_matrix)+1)]
    sM[0,0] = 0
    pX[0,:] = 1 #'X'
    pY[:,0] = 2 #'Y'
    pY[1,0] = pX[0,1] = 3 #'M'
    
    
    # Filling the matrices
    for i in range(1, len(s_query)+1):
        for j in range(1, len(s_reference_matrix)+1):
            probas_ref = s_reference_matrix[j-1]
            letter_query = s_query_one_hot[i-1]
            sub_score = ungapped_extension.UNGAPPED_SCORE_ALGORITHM[score_method](probas_ref, 
                                                                                  letter_query,
                                                                                  mismatch_score,
                                                                                  substitution)
            # First M
            M1 = sM[i-1, j-1] + sub_score
            X1 = sX[i-1, j-1] + sub_score
            Y1 = sY[i-1, j-1] + sub_score
            sM[i,j] = np.max([M1, X1, Y1])
            pM[i,j] = [3,1,2][np.argmax([M1, X1, Y1])]
            
            # Now X
            M2 = sM[i, j-1] + gap_bias + gap_penalty
            X2 = sX[i, j-1] + gap_penalty
            Y2 = sY[i, j-1] + gap_bias + gap_penalty
            sX[i,j] = np.max([M2, X2, Y2])
            pX[i,j] = [3,1,2][np.argmax([M2, X2, Y2])]
            
            # Finally Y
            M3 = sM[i-1, j] + gap_bias + gap_penalty
            X3 = sX[i-1, j] + gap_bias + gap_penalty
            Y3 = sX[i-1, j] + gap_penalty
            sY[i,j] = np.max([M3, X3, Y3])
            pY[i,j] = [3,1,2][np.argmax([M3, X3, Y3])]
    
    
    # Traceback
    # Starting from the bottom right...
    k = len(s_query)
    l = len(s_reference_matrix)
    
    max_score = np.max([sX[k,l], sY[k,l], sM[k,l]])
    
    if sX[k,l] == max_score:
        origin = pX[k,l]
    elif sY[k,l] == max_score:
        origin = pY[k,l]
    else :
        origin = pM[k,l]
    
    alignment_length = l
    string_query = ""
    string_reference = ""
    
    Z = "AXYM"
    path1 = [Z[int(origin)]]
    while (k > 0) | (l > 0):
        if origin == 2:
            origin = pY[(k, l)]
            k -= 1
        elif origin == 1:
            origin = pX[(k, l)]
            l -= 1
        elif origin == 3:
            origin = pM[(k, l)]
            l -= 1
            k -= 1
        if (k > 0) | (l > 0):
            path1.append(Z[int(origin)])
    
    # TODO optimize this and address case reverse = True
    path1 = path1[::-1]
    print(path1)
    string_query1 = ""
    string_reference1 = ""
    query_idx = 0
    ref_idx = 0
    
    
    for step in path1:
        print(string_query1)
        print(query_idx)
        print(ref_idx)
        if step == "X":
            string_reference1 += "ACGT"[np.argmax(s_reference_matrix[ref_idx])]
            string_query1 += "-"
            ref_idx += 1
        elif step == "Y":
            string_reference1 += "-"
            string_query1 += query[query_idx]
            query_idx += 1
        else:
            string_reference1 += "ACGT"[np.argmax(s_reference_matrix[l-1])]
            string_query1 += query[query_idx]
            ref_idx += 1
            query_idx += 1
    
    #If there are gaps at the end of the query, we remove them
    to_remove = 0
    end = len(string_query1)-to_remove
    start = len(string_query1)-to_remove-1
    while string_query1[start:end] == "-":
        to_remove += 1
        end = len(string_query1)-to_remove
        start = len(string_query1)-to_remove-1
    string_reference1 = string_reference1[0:(len(string_reference1) - to_remove)]
    string_query1 = string_query1[0:(len(string_query1) - to_remove)]
    ref_idx -= to_remove
    if to_remove > 0:
        max_score = max_score - gap_bias - gap_penalty*to_remove
    
    
    if reverse:        
        final_idx = ref_idx - pos_r
        string_reference1 = string_reference1[::-1]
        string_query1 = string_query1[::-1]
    else:
        final_idx = ref_idx + pos_r
    
    
    return(string_query1, string_reference1, final_idx, max_score)


def gapped_extension(query, reference, ungapped_dict, score_method,
                     substitution_dict, gap_penalty, gap_bias, mismatch_score):
    """
    compute a gapped alignment
    generates an alignment between the query and the reference sequences
    
    :param query: query sequence of letters (ACGT)
    :param reference: reference matrix file
    :param ungapped_dict: ungapped_extensions dict (positions and scores)
    :param score_method: method used to compute extension score
    :param substitution_dict: stores the cost of replacing one letter by another (typically from BLOSUM matrices)
    :param gap_bias: cost for opening a gap, if affine function
    :param gap_penalty: cost of a gap in the alignment
    
    :return: gapped_extensions dict (positions, scores and strings)
    """
    gapped_extensions = dict()
    
    for extension in ungapped_dict:
        # (tmp_pos_left_query, tmp_pos_right_query)] = [(tmp_pos_left_ref, tmp_pos_right_ref), tmp_score_left + tmp_score_right
        pos_left_q = extension[0]
        pos_right_q = extension[1]
        
        print(pos_left_q, pos_right_q)
        
        for pos_ref, score in ungapped_dict[extension]:
            # Extend to the left (only if there is something to extend)
            if pos_left_q > 0:
                (string_ql, string_rl, pos_l, score_l) = gapped_alignment(query, reference, pos_left_q, pos_ref[0], score_method, substitution_dict, mismatch_score, gap_penalty, reverse = True)
            else:
                string_ql = ""
                string_rl = ""
                pos_l = pos_ref[0]
                score_l = 0
            
            # Extend to the right  (only if there is something to extend)
            if pos_right_q < (len(query)-1):
                (string_qr, string_rr, pos_r, score_r) = gapped_alignment(query, reference, pos_right_q, pos_ref[1], score_method, substitution_dict, mismatch_score, gap_penalty)
            else:
                string_qr = ""
                string_rr = ""
                pos_r = pos_ref[1]
                score_l = 0
            
            new_score = score + score_l + score_r
            new_query_string = string_ql + query[pos_left_q:(pos_right_q+1)] + string_qr
            #TODO write a method to get the consensus sequence
            new_ref_string = string_rl + "".join(["ACGT"[np.argmax(s_reference_matrix[k])] for k in range(pos_ref[0],(pos_ref[1]+1))]) + string_rr
            gapped_extensions[(pos_l, pos_r)] = [new_query_string, new_ref_string, new_score]
    return gapped_extensions

# %%
'''
if __name__ == "main":

    def to_matrix(ref):
        mat = np.array([[0,0,0,0]])
        for i in ref:
            if i == "A":
                mat = np.append(mat, [[1, 0, 0, 0]], axis=0)
            elif i == "C":
                mat = np.append(mat, [[0, 1, 0, 0]], axis=0)
            elif i == "G":
                mat = np.append(mat, [[0, 0, 1, 0]], axis=0)
            else :
                mat = np.append(mat, [[0, 0, 0, 1]], axis=0)
        return(mat[1:])
    
    query = "AACTAATTTCCCGGGGATGAC"
    reference = "AAAACTAATTTCCCGTCGGAGGACTGC"
    reference_matrix = to_matrix(reference)
    
    gapped_alignment(query, reference, 15, 17, 'sum_proba_score', dict(), -1, gap_bias = 0, reverse = True)
    '''
# %%

query = "AAT"
ref = "AAA"

reference = np.zeros((3, 4))
reference[(0,0)] = 1
reference[(1,0)] = 1
reference[(2,0)] = 1


mismatch_score = 5

"""
expect: AAT VS AA-A ou AA-T VS AAA
"""

gapped_alignment(query, reference, -1, -1, 'sum_proba_score', substitution=0,
                     mismatch_score=mismatch_score, gap_penalty=-1, gap_bias = 0, reverse = False)