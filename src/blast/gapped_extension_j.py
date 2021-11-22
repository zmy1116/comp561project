import numpy as np
import pickle
from src.blast.nw_proba import nw_affine_two




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
    reference_matrix = pickle.load(open(reference, 'rb'))
    # reference_matrix = reference
    max_length = 3 * len(query)

    for extension in ungapped_dict:
        # (tmp_pos_left_query, tmp_pos_right_query)] = [(tmp_pos_left_ref, tmp_pos_right_ref), tmp_score_left + tmp_score_right
        pos_left_q = extension[0]
        pos_right_q = extension[1]

        do_left = True
        do_right = True

        if pos_left_q == 0:
            string_ql = ""
            string_rl = ""
            score_l = 0
            do_left = False

        if pos_right_q == len(query) - 1:
            string_qr = ""
            string_rr = ""
            score_r = 0
            do_right = False

        # print(pos_left_q, pos_right_q, do_left, do_right)


        for pos_ref, score in ungapped_dict[extension]:
            pos_l = pos_ref[0]
            pos_r = pos_ref[1]
            if do_right:
                pos_right_ref = pos_ref[1]
                if pos_right_ref < len(reference_matrix) - 1:
                    seq1 = query[pos_right_q + 1:]
                    seq2 = reference_matrix[pos_right_ref+1:min(len(reference_matrix),
                                            pos_right_ref+1+max_length)]
                    alignment, score_r = nw_affine_two(seq1, seq2, gap_bias,
                                                       gap_penalty, score_method,
                                                       substitution_dict,
                                                       offset=0, all_paths=False,
                                                       mismatch_score=mismatch_score)
                    string_qr = alignment[0]
                    string_rr = alignment[1]

                    # No need to consider the gaps at the end of the query

                    last_pos_r = len(string_qr) - 1
                    to_remove = 0
                    if  string_qr[last_pos_r] == "-": score_r -= gap_bias

                    while string_qr[last_pos_r] == "-":
                        to_remove += 1
                        last_pos_r -= 1
                        score_r -= gap_penalty

                    string_qr = string_qr[:last_pos_r+1]
                    string_rr = string_rr[: len(string_rr) - to_remove]

                    pos_r = pos_right_ref + len(string_rr) - string_rr.count("-")

                else:
                    # Align right part of the query with gaps
                    string_qr = query[pos_right_q + 1:]
                    string_rr = "-" * len(string_qr)
                    score_r = gap_bias + gap_penalty * len(string_qr)



            if do_left:
                pos_left_ref = pos_ref[0]
                if pos_left_ref > 0:
                    seq1 = query[:pos_left_q][::-1]
                    seq2 = reference_matrix[max(0, pos_left_ref - 1 - max_length):
                                            pos_left_ref][::-1]
                    alignment, score_l = nw_affine_two(seq1, seq2, gap_bias,
                                                       gap_penalty, score_method,
                                                       substitution_dict,
                                                       offset=0, all_paths=False,
                                                       mismatch_score=mismatch_score)


                    string_ql = alignment[0][::-1]
                    string_rl = alignment[1][::-1]

                    # No need to consider the gaps at the beginning of the query

                    last_pos_l = 0
                    if  string_ql[last_pos_l] == "-": score_l -= gap_bias

                    while string_ql[last_pos_l] == "-":
                        last_pos_l += 1
                        score_l -= gap_penalty

                    string_ql = string_ql[last_pos_l:]
                    string_rl = string_rl[last_pos_l:]

                    pos_l = pos_left_ref - len(string_rl) + string_rl.count("-")

                else:
                    # Align left part of the query with gaps
                    string_ql = query[0 : pos_left_q]
                    string_rl = "-" * len(string_ql)
                    score_l = gap_bias + gap_penalty * len(string_ql)

            new_score = score + score_l + score_r
            new_query_string = string_ql + query[pos_left_q:(pos_right_q+1)] + string_qr
            new_ref_string = string_rl + "".join(["ACGT"[np.argmax(reference_matrix[k])] for k in range(pos_ref[0],(pos_ref[1]+1))]) + string_rr
            gapped_extensions[(pos_l, pos_r)] = [new_query_string, new_ref_string, new_score]

    return gapped_extensions

# %%

"""
query = "AATA"
ref = "AAAA"

reference = np.zeros((4, 4))
reference[(0,0)] = 1
reference[(1,0)] = 1
reference[(2,0)] = 1
reference[(3,0)] = 1

mismatch_score = 5

"""
#expect: AAT VS AA-A ou AA-T VS AAA
"""

ungapped_dict=dict()
ungapped_dict[(3,3)] = [[(3,3), 0]]

score_method = "sum_proba_score"

gapped_extension(query, reference, ungapped_dict, score_method,
                 substitution_dict=dict(), gap_penalty=0,
                 gap_bias=0, mismatch_score=mismatch_score)

"""

# pos_right_ref + len(string_rr) - string_rr.count("-") - pos_left_ref +  len(string_rl) - string_rr.count("-")