import pickle
import numpy as np


def generate_queries(reference_matrix, query_length, type='all_positions', num=1, with_substitution=False,
                     characters="ACGT"):
    """
    generate testing sequences according to ref seq probabilities
    :param reference_matrix: matrix storing the proba of each letter at each
    position in ref sequence
    :param query_length: desired length of the query sequence
    :param type: generation type (all_positions: for every position generate queries, random: choose a random position)
    :param num: number of queries generated for each position
    :param with_substitution: whether perform extra substitutions
    :param characters: ordered list of characters, same order as in reference matrix
    :output: a query sequence generated according to the probabilities of the
    reference sequence, no extra substitution or indel yet
    """
    assert type in ['all_positions', 'random']

    ref_length, nb_chars = reference_matrix.shape
    # define the choices of query positions
    query_position_ranges = np.arange(0, ref_length - query_length)

    # select list of labels: all positions or just a random position
    if type == 'all_position':
        # iterate through all positions
        query_positions = query_position_ranges
    else:
        # choose a random position
        query_positions = [np.random.choice(query_position_ranges)]

    # for each query position, generate #num queries
    queries_data = {}
    for label in query_positions:
        queries = []
        for _ in range(num):
            query = ''
            for k in range(query_length):
                probas = reference_matrix[k + label]
                char_index = np.random.choice(np.arange(nb_chars), p=probas)
                query += characters[char_index]

            # TODO ADD SUBSTITUTION AND GAP
            # if with_substitution:
            #     query = add_substitutions(query, substitutions_probas=0, non_substitution_proba=0.9,
            #                               characters=characters)
            queries.append(query)
        queries_data[label] = queries

    return queries_data


def add_substitutions(query, substitutions_probas=0, non_substitution_proba=0.9, characters="ACGT"):
    """
    Starting from a query sequence, add substitutions to it
    :param query: query sequence to modify
    :param substitutions_probas: dictionary containing the probabilities of
    switching from one character to another, includes the proba of not changing
    in substitutions_probas[(c1, c1)], c1 being the index of corresponding character
    in characters
    :param non_substitution_proba: if no substitution dict is given, proba of not
    having a substitution
    :param characters: ordered list of characters, same order as in reference matrix
    :output: modified query with random substitutions in it

    """
    nb_chars = len(characters)

    if substitutions_probas == 0:
        # equiprobability for all substitutions
        substitutions_probas = dict()
        for c1 in characters:
            for c2 in characters:
                if c2 == c1:
                    substitutions_probas[(c1, c2)] = non_substitution_proba
                else:
                    substitutions_probas[(c1, c2)] = (1 - non_substitution_proba) / (len(characters) - 1)
    new_query = ""
    for i, c in enumerate(query):
        probas = [substitutions_probas[(c, c2)] for c2 in characters]
        char_index = np.random.choice(np.arange(nb_chars), p=probas)
        new_query = new_query + characters[char_index]

    return new_query
