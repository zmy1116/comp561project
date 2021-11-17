import numpy as np

def sequence_one_hot(seq, letters='ACGT'):
    """
    generate one hot encoding of input sequence
    :param seq: raw sequence
    :param letters:  letters in order
    :return: one hot encodoing matrix
    """
    letter_map = {
        letters[0]: [1, 0, 0, 0],
        letters[1]: [0, 1, 0, 0],
        letters[2]: [0, 0, 1, 0],
        letters[3]: [0, 0, 0, 1]
    }
    seq = [letter_map[x] for x in seq]
    seq = np.array(seq)
    return seq


def seq2num(seq, letters='ACGT'):
    """
    generate numerical representation of input sequence
    :param seq: raw sequence
    :param letters:  letters in order
    :return: numerical array
    """
    letter_map = {
        letters[0]: 0,
        letters[1]: 1,
        letters[2]: 2,
        letters[3]: 3
    }
    seq = [letter_map[x] for x in seq]
    seq = np.array(seq)
    return seq


def generate_sequence(possible_letters_dict, fixed_length,
                      saved_strings, current_string, characters="ACGT"):
    """
    Recursively find all possible strings of fixed length with letter at rank
    k taken among letters in possible_letters_dict[k]
    :param possible_letters_dict: possible letters for each position
    :fixed_length: length of strings to generate
    :param saved_strings: set of strings already computed
    :param current_string: string being completed
    :param characters: characters to use, order is the same as indices in
    possible_letters_dict
    :return: nothing but fills saved_strings
    """

    if len(current_string) == fixed_length:
        saved_strings.add(current_string)
    else:
        pos = len(current_string)
        possible_letters = possible_letters_dict[pos]
        for idx in possible_letters:
            new_current_string = current_string + characters[idx]
            generate_sequence2(possible_letters_dict, fixed_length, saved_strings,
                              new_current_string, characters="ACGT")


def generate_sequence3(possible_letters_dict, fixed_length, reference_matrix,
                       saved_scores, saved_strings, current_score,
                       current_string, characters="ACGT"):
    """
    Recursively find all possible strings of fixed length with letter at rank
    k taken among letters in possible_letters_dict[k]
    :param possible_letters_dict: possible letters for each position
    :fixed_length: length of strings to generate
    :param saved_strings: set of strings already computed
    :param current_string: string being completed
    :param characters: characters to use, order is the same as indices in
    possible_letters_dict
    :return: nothing but fills saved_strings
    """

    if len(current_string) == fixed_length:
        saved_strings.add(current_string)
        saved_scores[current_string] = current_score

    else:
        pos = len(current_string)
        possible_letters = possible_letters_dict[pos]
        for idx in possible_letters:
            new_current_string = current_string + characters[idx]
            new_current_score = current_score + np.log(reference_matrix[pos][idx])
            generate_sequence3(possible_letters_dict, fixed_length,
                               reference_matrix, saved_scores, saved_strings,
                               new_current_score, new_current_string,
                               characters="ACGT")

# %%
"""
possible_letters_dict = above_threshold
fixed_length = 11
reference_matrix = reference_matrix[0:11]
current_string = ("", 0, "")
saved_strings = set()

generate_sequence(possible_letters_dict, fixed_length, R,
                      saved_strings, current_string, characters="ACGT")
"""