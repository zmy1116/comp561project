import numpy as np
import pickle
from src.utils.utils import seq2num
from src.utils.registry import Registry
from src.utils.utils import generate_sequence

HASHTABLE_SEEDING_ALGORITHM = Registry()
HASHTABLE_MATCHING_ALGORITHM = Registry()


@HASHTABLE_SEEDING_ALGORITHM.register('consensus_seed_seq')
def consensus_seq_method(reference_matrix, k=11):
    """
    method 1
    take the most likely sequence out of reference matrix, and directly building seed table

    :param reference_matrix: reference probability matrix
    :param k:  seed size
    :return:  seed table
    """
    consensus_seq = np.argmax(reference_matrix, axis=1)

    table = {}

    for idx in range(len(consensus_seq) - k + 1):

        seed = tuple(consensus_seq[idx:idx + k])
        if seed in table:
            table[seed].append(idx)
        else:
            table[seed] = [idx]

    return table



@HASHTABLE_SEEDING_ALGORITHM.register('individual_nt_threshold')
def individual_nt_threshold(reference_matrix, k=11, threshold=0.15,
                            characters="ACGT"):
    """
    method 1
    take the most likely sequence out of reference matrix, and directly building seed table

    :param reference_matrix: reference probability matrix
    :param k:  seed size
    :param threshold:  min proba to take into account a nucleotide
    :return:  seed table
    """
    N = reference_matrix.shape[0]

    table = {}

    above_threshold = dict()

    for idx in range(N):
        line_above_threshold = list(np.where(reference_matrix[idx] >= threshold)[0])
        above_threshold[idx] = line_above_threshold

    current_seeds = set()
    generate_sequence(above_threshold, k, current_seeds,
                      "", characters="ACGT")

    for seed in current_seeds: table[seed] = [0]

    for idx in range(1, N-k+1):
        new_current_seeds = {seed[1:] for seed in current_seeds}
        current_seeds = set()
        for char_idx in above_threshold[idx + k - 1]:
            current_seeds = current_seeds.union({s + characters[char_idx]
                                                 for s in new_current_seeds})
        for seed in current_seeds:
            if seed in table:
                table[seed].append(idx)
            else:
                table[seed] = [idx]

    return table


def hashtable_generation(reference_matrix_file, method, k, output_path=""):
    """
    generate and store seed table for blast
    :param reference_matrix_file: reference matrix file
    :param method: table generation method
    :param k: seed size
    :param output_path: output table storage path
    :return: hashtable
    """
    assert method in HASHTABLE_SEEDING_ALGORITHM
    reference_matrix = pickle.load(open(reference_matrix_file, 'rb'))
    table = HASHTABLE_SEEDING_ALGORITHM[method](reference_matrix, k)
    table_data = {
        'method': method,
        'hashtable': table,
        'k': k
    }
    if output_path != "": pickle.dump(table_data, open(output_path, 'wb'))
    return table_data


def consensus_seq_match(table_data, query):
    """

    :param table_data: hashtable data
    :param query: query sequence of letters (ACGT)
    :return: list of dictionary with key query_start_idx,  reference_start_idx, match_score
    """
    hashtable = table_data['hashtable']
    k = table_data['k']
    results = []
    query = seq2num(query)
    for idx in range(len(query) - k + 1):
        query_seed = tuple(query[idx:idx+k])
        matches = hashtable.get(query_seed, [])
        if len(matches) > 0:
            results.append(
                {
                    'query_idx': idx,
                    'db_indices': matches
                }
            )
    return results
