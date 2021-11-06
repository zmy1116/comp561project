import pickle
import numpy as np


def generate_vanilla_query(reference_matrix, query_length, characters="ACGT"):
    """
    start from random position in ref sequence and generate query seq according
    to ref seq probabilities
    :param reference_matrix_file: matrix storing the proba of each letter at each
    position in ref sequence
    :param length: desired length of the query sequence
    :param characters: ordered list of characters, same order as in reference matrix
    :output: a query sequence generated according to the probabilities of the
    reference sequence, no extra substitution or indel yet
    """
    # reference_matrix = pickle.load(open(reference_matrix_file, 'rb'))
    ref_length, nb_chars = reference_matrix.shape
    query_pos = np.random.choice(np.arange(0, ref_length - query_length))
    print(query_pos)
    query = ""
    for k in range(query_length):
        probas = reference_matrix[query_pos]
        char_index = np.random.choice(np.arange(nb_chars), p=probas)
        query += characters[char_index]
        query_pos += 1
    return query


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
                if c2 == c1: substitutions_probas[(c1, c2)]=non_substitution_proba
                else:
                    substitutions_probas[(c1, c2)] = (1-non_substitution_proba) / (len(characters)-1)
    new_query = ""
    for i, c in enumerate(query):
        probas = [substitutions_probas[(c, c2)] for c2 in characters]
        char_index = np.random.choice(np.arange(nb_chars), p=probas)
        new_query = new_query + characters[char_index]

    return new_query







