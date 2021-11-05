from src.utils.utils import sequence_one_hot
from src.utils.registry import Registry

UNGAPPED_SCORE_ALGORITHM = Registry()

@UNGAPPED_SCORE_ALGORITHM.register('sum_proba_score')
def sum_proba_score(ref_letter_prob, query_letter_one_hot, substitution=dict()):
    """
    compute score for
    :param ref_letter_prob: 1D array giving the proba for each letter in ref
    :param query_letter_one_hot: one hot encoding of the letter in query
    :param substitution: dict storing cost of replacing one letter with another
    :output: score of matching query_letter_one_hot at position corresponding to ref_letter_prob
    """
    score = 0
    for i, p in enumerate(ref_letter_prob):
        score += p * query_letter_one_hot[i] - p * (query_letter_one_hot[i] != 1)
    return score


def ungapped_extension(query, matches_dict, reference_matrix_file, k, delta,
                       score_method, substitution=dict()):
    """
    compute ungapped extension
    generates ungapped extended matches and corresponding HSP scores
    :param query: query sequence of letters (ACGT)
    :param reference_matrix_file: reference matrix file
    :matches: matches between query seq and seeds, see consensus_seq_match for format
    :param k: seed size
    :param delta: max allowed score drop before extension is stopped, positive value
    :param score_method: method used to compute extension score
    :param substitution_dict: stores the cost of replacing one letter by another (typically from BLOSUM matrices)
    :return: ungapped_extensions dict (positions and scores)
    """
    assert score_method in UNGAPPED_SCORE_ALGORITHM

    query_one_hot = sequence_one_hot(query)
    reference_matrix = pickle.load(open(reference_matrix_file, 'rb'))

    ungapped_extensions = dict()

    for matches in matches_dict:
        query_left = matches['query_idx']
        query_right = query_left + k - 1

        for match in matches['db_indices']:
            tmp_score_left = 0
            current_score_left = 0
            tmp_pos_left_ref = match
            current_pos_left_ref = match
            tmp_pos_left_query = query_left
            current_pos_left_query = query_left

            tmp_score_right = 0
            current_score_right = 0
            tmp_pos_right_ref = match + k - 1
            current_pos_right_ref = match + k - 1
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
                current_score_left += UNGAPPED_SCORE_ALGORITHM[score_method](probas_ref, letter_query, substitution)

                if current_score_left > tmp_score_left:
                    tmp_score_left = current_score_left
                    tmp_pos_left_ref = current_pos_left_ref
                    tmp_pos_left_query = current_pos_left_query

                diff_left = current_score_left - tmp_score_left

            # Right
            diff_right = current_score_right - tmp_score_right
            while (diff_right > -delta) and (current_pos_right_ref > 0) and (current_pos_right_query > 0):
                current_pos_right_ref += 1
                current_pos_right_query += 1

                probas_ref = reference_matrix[current_pos_right_ref]
                letter_query = query_one_hot[current_pos_right_query]

                current_score_right += UNGAPPED_SCORE_ALGORITHM[score_method](probas_ref, letter_query, substitution)

                if current_score_right > tmp_score_right:
                    tmp_score_right = current_score_right
                    tmp_pos_right_ref = current_pos_right_ref
                    tmp_pos_right_query = current_pos_right_query

                diff_right = current_score_right - tmp_score_right

            # Store results
            ungapped_extensions[(tmp_pos_left_query, tmp_pos_right_query)] = [(tmp_pos_left_ref, tmp_pos_right_ref), tmp_score_left + tmp_score_right]

    return ungapped_extensions