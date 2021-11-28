import numpy as np
import pickle
from src.blast.nw_proba import nw_affine_two


def clean_end_gaps(string_query, string_ref, score, gap_bias, gap_penalty):
    """
    remove end gaps and modify scores if necessary
    clean end side gaps
    :param string_query:
    :param string_ref:
    :param score:
    :param gap_bias: cost for opening a gap, if affine function
    :param gap_penalty: cost of a gap in the alignment
    :return: query string, reference string, score
    """

    if string_query[-1] == '-':
        score -= gap_bias

    while string_query[-1] == '-':
        string_query = string_query[:-1]
        string_ref = string_ref[:-1]
        score -= gap_penalty

    return string_query, string_ref, score


def gapped_extension_one_side(seq_query, seq_ref, score_method,
                              substitution_dict, gap_penalty, gap_bias, mismatch_score, side):
    if side == 'left':
        seq_query = seq_query[::-1]
        seq_ref = seq_ref[::-1]

    if len(seq_ref) == 0:
        string_query = seq_query
        string_ref = '-' * len(string_query)
        score = gap_bias + gap_penalty * len(string_query)
        aligned = False
    else:
        [string_query, string_ref], score = nw_affine_two(seq_query, seq_ref, gap_bias,
                                                          gap_penalty, score_method,
                                                          substitution_dict,
                                                          offset=0, all_paths=False,
                                                          mismatch_score=mismatch_score)
        try:
            string_query, string_ref, score = clean_end_gaps(string_query, string_ref, score, gap_bias,
                                                             gap_penalty)
        except:
            print(string_query, string_ref)
            raise
        aligned = True

    if side == 'left':
        string_query = string_query[::-1]
        string_ref = string_ref[::-1]

    return string_query, string_ref, score, aligned


def gapped_extension(query, reference, ungapped_dict, score_method,
                     substitution_dict, gap_penalty, gap_bias, mismatch_score, ref_max_length_factor):
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
    :param ref_max_length_factor:  do gapped extention of length FACTOR x ref_max_length_factor
    :return: gapped_extensions dict (positions, scores and strings)
    """

    reference_matrix = pickle.load(open(reference, 'rb'))
    # reference_matrix = reference
    max_length = ref_max_length_factor * len(query)

    gapped_extensions = []

    for extension in ungapped_dict:

        ungapped_entry = extension['ungapped_extension_result']

        # query left and right matching indices
        pos_left_q = ungapped_entry['query_left_idx']
        pos_right_q = ungapped_entry['query_right_idx']
        do_left = pos_left_q != 0
        do_right = pos_right_q != (len(query) - 1)

        # reference left and right matching indices
        pos_l = ungapped_entry['ref_left_idx']
        pos_r = ungapped_entry['ref_right_idx']

        # right side
        string_query, string_ref, score = '', '', 0
        ref_aligned_indices = -1, -1
        if do_right:
            # get reference and query sequences to be aligned
            ref_boundary = pos_r + 1, pos_r + 1 + max_length
            seq_ref = reference_matrix[ref_boundary[0]: ref_boundary[1]]
            seq_query = query[pos_right_q + 1:]
            string_query, string_ref, score, aligned = gapped_extension_one_side(seq_query, seq_ref, score_method,
                                                                                 substitution_dict, gap_penalty,
                                                                                 gap_bias, mismatch_score, 'right')
            if aligned:
                ref_aligned_indices = pos_r + 1, pos_r + len(string_ref) - string_ref.count("-")
        right_alignment_result = {
            'ref_left_idx': ref_aligned_indices[0],
            'ref_right_idx': ref_aligned_indices[-1],
            'score': score,
            'query_string': string_query,
            'ref_string': string_ref
        }

        # left side
        string_query, string_ref, score = '', '', 0
        ref_aligned_indices = -1, -1
        if do_left:
            ref_boundary = max(0, pos_l - 1 - max_length), pos_l
            seq_ref = reference_matrix[ref_boundary[0]:ref_boundary[1]]
            seq_query = query[:pos_left_q]

            string_query, string_ref, score, aligned = gapped_extension_one_side(seq_query, seq_ref, score_method,
                                                                                 substitution_dict, gap_penalty,
                                                                                 gap_bias, mismatch_score, 'left')
            if aligned:
                ref_aligned_indices = pos_l - (len(string_ref) - string_ref.count("-")), pos_l - 1

        left_alignment_result = {
            'ref_left_idx': ref_aligned_indices[0],
            'ref_right_idx': ref_aligned_indices[-1],
            'score': score,
            'query_string': string_query,
            'ref_string': string_ref
        }

        # aggregate
        if left_alignment_result['ref_left_idx'] == -1:
            ref_left_idx = pos_l
        else:
            ref_left_idx = left_alignment_result['ref_left_idx']
        ref_right_idx = max(pos_r, right_alignment_result['ref_right_idx'])
        gapped_extension_result = {
            'ref_left_idx': ref_left_idx,
            'ref_right_idx': ref_right_idx,
            'score': left_alignment_result['score'] + right_alignment_result['score'],
            'left_alignment_result': left_alignment_result,
            'right_alignment_result': right_alignment_result,
        }

        extension['gapped_extension_result'] = gapped_extension_result
        final_result = {
            'ref_left_idx': gapped_extension_result['ref_left_idx'],
            'ref_right_idx': gapped_extension_result['ref_right_idx'],
            'score': gapped_extension_result['score'] + extension['ungapped_extension_result']['score'] +
                     extension['seed_matching_result']['score']
        }
        extension['final_result'] = final_result

        gapped_extensions.append(extension)

    return gapped_extensions
