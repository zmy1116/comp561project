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


def generate_sequence(possible_letters_dict, fixed_length, saved_strings,
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
    if len(current_string) == fixed_length: saved_strings.add(current_string)
    else:
        pos = len(current_string)
        possible_letters = possible_letters_dict[pos]
        for idx in possible_letters:
            new_current_string = current_string + characters[idx]
            generate_sequence(possible_letters_dict, fixed_length, saved_strings,
                              new_current_string, characters="ACGT")

