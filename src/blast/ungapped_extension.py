import pickle
from src.utils.utils import sequence_one_hot
from src.blast.nt_scoring_function import NT_SCORE_ALGORITHM


def ungapped_extension(query, matches_dict, reference_matrix_file, k, delta,
                       score_method, mismatch_score=1, substitution=dict()):
    """
    compute ungapped extension
    generates ungapped extended matches and corresponding HSP scores
    :param query: query sequence of letters (ACGT)
    :param reference_matrix_file: reference matrix file
    :matches: matches between query seq and seeds, see consensus_seq_match for format
    :param k: seed size
    :param delta: max allowed score drop before extension is stopped, positive value
    :param score_method: method used to compute extension score
    :param substitution_dict: stores the cost of replacing one letter by another
    :return: ungapped_extensions dict (positions and scores)
    """
    assert score_method in NT_SCORE_ALGORITHM

    query_one_hot = sequence_one_hot(query)
    reference_matrix = pickle.load(open(reference_matrix_file, 'rb'))
    ref_length, _ = reference_matrix.shape

    ungapped_extensions = []

    for match in matches_dict:
        query_idx = match['seed_matching_result']['query_idx']
        match_idx = match['seed_matching_result']['ref_idx']
        query_left = query_idx
        query_right = query_left + k - 1

        tmp_score_left = 0
        current_score_left = 0
        tmp_pos_left_ref = match_idx
        current_pos_left_ref = match_idx
        tmp_pos_left_query = query_left
        current_pos_left_query = query_left

        tmp_score_right = 0
        current_score_right = 0
        tmp_pos_right_ref = match_idx + k - 1
        current_pos_right_ref = match_idx + k - 1
        tmp_pos_right_query = query_right
        current_pos_right_query = query_right

        # Separately treat left and right extensions

        # Left
        diff_left = current_score_left - tmp_score_left
        while (diff_left > -delta) and (current_pos_left_ref > 0) and (current_pos_left_query > 0):
            current_pos_left_ref -= 1
            current_pos_left_query -= 1

            probas_ref = reference_matrix[current_pos_left_ref]
            letter_query = query_one_hot[current_pos_left_query]
            current_score_left += NT_SCORE_ALGORITHM[score_method](probas_ref, letter_query, mismatch_score,
                                                                   substitution)

            if current_score_left > tmp_score_left:
                tmp_score_left = current_score_left
                tmp_pos_left_ref = current_pos_left_ref
                tmp_pos_left_query = current_pos_left_query

            diff_left = current_score_left - tmp_score_left

        # Right
        diff_right = current_score_right - tmp_score_right
        while (diff_right > -delta) and (current_pos_right_ref < ref_length - 1) and (
                current_pos_right_query < len(query) - 1):
            current_pos_right_ref += 1
            current_pos_right_query += 1

            probas_ref = reference_matrix[current_pos_right_ref]
            letter_query = query_one_hot[current_pos_right_query]

            current_score_right += NT_SCORE_ALGORITHM[score_method](probas_ref, letter_query, mismatch_score,
                                                                    substitution)

            if current_score_right > tmp_score_right:
                tmp_score_right = current_score_right
                tmp_pos_right_ref = current_pos_right_ref
                tmp_pos_right_query = current_pos_right_query

            diff_right = current_score_right - tmp_score_right

        extension_result = {
            'query_left_idx': tmp_pos_left_query,
            'query_right_idx': tmp_pos_right_query,
            'ref_left_idx': tmp_pos_left_ref,
            'ref_right_idx': tmp_pos_right_ref,
            'score': tmp_score_left + tmp_score_right
        }
        match['ungapped_extension_result'] = extension_result
        ungapped_extensions.append(match)

    return ungapped_extensions
