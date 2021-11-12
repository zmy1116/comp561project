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
    reference_matrix = pickle.load(open(reference, 'rb'))
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
        s_query = query[(pos_q+1):]
        end = np.min([pos_r+max_length, len(reference_matrix)])
        s_reference_matrix = reference_matrix[(pos_r+1):end]
        
    s_query_one_hot = sequence_one_hot(s_query)
    sM = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    sX = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    sY = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pM = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pX = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    pY = np.zeros((len(s_query)+1, len(s_reference_matrix)+1))
    
    print(pY.shape)
    
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
            sub_score = ungapped_extension.UNGAPPED_SCORE_ALGORITHM[score_method](probas_ref, letter_query, mismatch_score, substitution)
            M1 = sM[i-1, j-1] + sub_score
            X1 = sX[i-1, j-1] + sub_score
            Y1 = sY[i-1, j-1] + sub_score
            sM[i,j] = np.max([M1, X1, Y1])
            pM[i,j] = [3,1,2][np.argmax([M1, X1, Y1])]
            M2 = sM[i-1, j] + gap_bias + gap_penalty
            X2 = sX[i-1, j] + gap_penalty
            sX[i,j] = np.max([M2, X2])
            pX[i,j] = [3,1][np.argmax([M2, X2])]
            M3 = sM[i, j-1] + gap_bias + gap_penalty
            Y3 = sX[i, j-1] + gap_penalty
            sY[i,j] = np.max([M3, Y3])
            pY[i,j] = [3,2][np.argmax([M3, Y3])]
    '''
            score_l = scores_mat[i-1,j] + gap_penalty # score from the left
            score_a = scores_mat[i,j-1] + gap_penalty # score from above
            
            # score from the diagonal
            probas_ref = s_reference_matrix[j]
            letter_query = s_query_one_hot[i]
            score_d = scores_mat[i-1,j-1] + UNGAPPED_SCORE_ALGORITHM[score_method](probas_ref, letter_query, substitution)
            
            scores_mat[i,j] = np.max([score_l, score_a, score_d])
            trace_mat[i,j] = np.argmax([score_l, score_a, score_d])
    '''
    
    # Traceback
    # Starting from the bottom right...
    k = len(s_query)
    l = len(s_reference_matrix)
    
    #If we have gaps at the end of the query, we ignore them and take the best score after the gaps
    x = l
    while pX[k,x] == 1:
        x =- 1
    max_score = np.max([sX[k,x], sY[k,l], sM[k,l]])
    if sX[k,x] == max_score:
        l = x
        origin = pX[k,l]
    elif sY[k,l] == max_score:
        origin = pY[k,l]
    else :
        origin = pM[k,l]
    
    alignment_length = l
    string_query = ""
    string_reference = ""
    
    while (k > 0) | (l > 0):
        if origin == 1: # From the left
            string_query = "-" + string_query
            letter = "ACGT"[np.argmax(s_reference_matrix[l-1])]
            string_reference = letter + string_reference
            l -= 1
            origin = pX[k,l]
        elif origin == 2: # From above
            string_query = s_query[k-1] + string_query
            string_reference = "-" + string_reference
            k -= 1
            origin = pY[k,l]
        else : # From the diagonal
            string_query = s_query[k-1] + string_query
            letter = "ACGT"[np.argmax(s_reference_matrix[l-1])]
            string_reference = letter + string_reference
            l -= 1
            k -= 1
            origin = pM[k,l]
        
    '''
    # We don't care about the gaps at the end of the query, so we loop until we have none
    while trace_mat[k,l] == 0:
        l -= 1
    
    final_score = scores_mat[k,l]
    string_query = ""
    string_reference = ""
    # Until we reach the top left, trace back and record the alignment
    while (k > 0) | (l > 0):
        # If we come from the left, add a gap to the reference
        if trace_mat[k,l] == 0:
            string_query = s_query[k] + string_query
            string_reference = "-" + string_reference
            l -= 1
        # If we come from above, add a gap to the query
        elif trace_mat[k,l] == 1:
            string_query = "-" + string_query
            letter = "ACGT"[np.argmax(s_reference_matrix[l])]
            string_reference = letter + string_reference
            k -= 1
        # If we come from the diagonal, add a the letters to the reference
        else :
            string_query = s_query[k] + string_query
            letter = "ACGT"[np.argmax(s_reference_matrix[l])]
            string_reference = letter + string_reference
            l -= 1
            k -= 1
    '''
    # Reverting if necessary
    if reverse:
        final_ref_position = pos_r - alignment_length
        string_query = string_query[::-1]
        string_reference = string_reference[::-1]
    else:
        final_ref_position = pos_r + alignment_length    
    return(string_query, string_reference, final_ref_position, final_score)


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
    #
    COUNT = 0
    #
    
    for extension in ungapped_dict:
        # (tmp_pos_left_query, tmp_pos_right_query)] = [(tmp_pos_left_ref, tmp_pos_right_ref), tmp_score_left + tmp_score_right
        pos_left_q = extension[0]
        pos_right_q = extension[1]
        
        print(pos_left_q, pos_right_q)
        
        for pos_ref, score in ungapped_dict[extension]:
            # Extend to the left
            (string_ql, string_rl, pos_l, score_l) = gapped_alignment(query, reference, pos_left_q, pos_ref[0], score_method, substitution_dict, mismatch_score, gap_penalty, reverse = True)
                
            # Extend to the right
            (string_qr, string_rr, pos_r, score_r) = gapped_alignment(query, reference, pos_right_q, pos_ref[1], score_method, substitution_dict, mismatch_score, gap_penalty)
                
            new_score = score + score_l + score_r
            new_query_string = string_ql + query[pos_left_q:(pos_right_q+1)] + string_qr
            #TODO write a method to get the consensus sequence
            new_ref_string = string_rl + "".join(["ACGT"[np.argmax(s_reference_matrix[k])] for k in range(pos_ref[0],(pos_ref[1]+1))]) + string_rr
            gapped_extensions[(pos_l, pos_r)] = [new_query_string, new_ref_string, new_score]
    return gapped_extensions

# %%

if __name__ == "main":
    '''
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
    '''
    gapped_alignment(query, reference, 15, 17, 'sum_proba_score', dict(), -1, gap_bias = 0, reverse = True)