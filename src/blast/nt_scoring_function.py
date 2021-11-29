"""
All nt scoring functions
"""


from src.utils.registry import Registry
import numpy as np

NT_SCORE_ALGORITHM = Registry()


@NT_SCORE_ALGORITHM.register('sum_proba_score')
def sum_proba_score(ref_letter_prob, query_letter_one_hot, mismatch_score=1,
                    substitution=dict()):
    """
    compute score for
    :param ref_letter_prob: 1D array giving the proba for each letter in ref
    :param query_letter_one_hot: one hot encoding of the letter in query
    :param substitution: dict storing cost of replacing one letter with another
    :output: score of matching query_letter_one_hot at position corresponding to ref_letter_prob
    """
    score = np.sum(
        ref_letter_prob * query_letter_one_hot - mismatch_score * ref_letter_prob * (1 - query_letter_one_hot))

    return score


@NT_SCORE_ALGORITHM.register('sum_proba_score_correct0')
def sum_proba_score_correct0(ref_letter_prob, query_letter_one_hot, mismatch_score=1,
                             substitution=dict()):
    """
    a modified scoring function, compute the probability score, when the letter matches the maximum probability letter,
    score = max(0, score)
    :param ref_letter_prob: 1D array giving the proba for each letter in ref
    :param query_letter_one_hot: one hot encoding of the letter in query
    :param substitution: dict storing cost of replacing one letter with another
    :output: score of matching query_letter_one_hot at position corresponding to ref_letter_prob
    """
    score = np.sum(
        ref_letter_prob * query_letter_one_hot - mismatch_score * ref_letter_prob * (1 - query_letter_one_hot))

    if np.argmax(ref_letter_prob) == np.argmax(query_letter_one_hot):
        score = max(score, 0)

    return score
