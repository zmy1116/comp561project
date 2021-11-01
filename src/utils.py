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