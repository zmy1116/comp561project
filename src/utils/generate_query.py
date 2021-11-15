import pickle
import numpy as np


def generate_queries(reference_matrix, query_length, query_positions, num=1,
                     with_substitution=False, with_indel=False, characters="ACGT"):
    """
    generate testing sequences according to ref seq probabilities
    :param reference_matrix: matrix storing the proba of each letter at each
    position in ref sequence
    :param query_length: desired length of the query sequence
    :param query_positions: starting position of the query, if it is a list => it's the list of query positions, if it's 0 => a randomly generated position, other number => all positions
    :param num: number of queries generated for each position
    :param with_substitution: whether perform extra substitutions
    :param characters: ordered list of characters, same order as in reference matrix
    :output: list of dictionary {label position, query string}, no extra substitution or indel yet
    """
    ref_length, nb_chars = reference_matrix.shape
    # define the choices of query positions
    query_position_ranges = np.arange(0, ref_length - query_length)

    # select list of labels: all positions or just a random position
    if type(query_positions) != list:

        if query_positions == 0:
            query_positions = [np.random.choice(query_position_ranges)]
        else:
            query_positions = list(query_position_ranges)

    # for each query position, generate #num queries
    queries_data = []
    for label in query_positions:
        for _ in range(num):
            query = ''
            for k in range(query_length):
                probas = reference_matrix[k + label]
                char_index = np.random.choice(np.arange(nb_chars), p=probas)
                query += characters[char_index]

            if with_substitution:
                query = add_substitutions(query, substitutions_probas=substitutions_probas,
                non_substitution_proba=0.9, characters=characters)
            print(s1)
            print(s1==query)
            # TODO ADD INDEL
            # if with_indel:
                # query = add_indel()

            queries_data.append(
                {
                    'pos': label,
                    'query': query
                }
            )
    return queries_data


# %% Substitutions probabilities
# Transitions are more likely than transversions

nucleotides = "ACGT"
non_substitution_proba = 0.9
substitutions_probas = dict()
purines = {"A", "G"}
pyrimidines = {"T", "C"}

for nt1 in nucleotides:
    for nt2 in nucleotides:
        if nt1 == nt2: substitutions_probas[(nt1, nt2)] = non_substitution_proba
        elif (nt1 in purines and nt2 in purines) or (nt1 in pyrimidines and nt2 in pyrimidines):
            substitutions_probas[(nt1, nt2)] = (1 - non_substitution_proba) / 2
        else: substitutions_probas[(nt1, nt2)] = (1 - non_substitution_proba) / 4



# %%
def add_substitutions(query, substitutions_probas=substitutions_probas, non_substitution_proba=0.9, characters="ACGT"):
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
