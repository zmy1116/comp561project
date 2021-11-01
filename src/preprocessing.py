import numpy as np
import pickle
import argparse
from src.utils import sequence_one_hot


def generate_reference_seq_matrix(seq_file, prob_file, output_path, letters='ACGT'):
    """
    generate reference sequence matrix size [#NT x 4]
    :param seq_file: file of most likely sequence
    :param prob_file: probability of most likely sequence
    :param letters: fix vocabulary in order
    :param output_path: saved data path
    :return:
    """
    # load sequence to one hot encoding
    seq = open(seq_file).readlines()[0]
    seq = sequence_one_hot(seq, letters)

    probs = open(prob_file).readlines()[0]
    probs = [eval(x) for x in probs.split(' ')[:-1]]
    probs = np.array(probs)

    rest_probs = (1 - probs) / 3
    seq_matrix = probs.reshape(-1, 1) * seq + rest_probs.reshape(-1, 1) * (1 - seq)
    pickle.dump(seq_matrix, open(output_path, 'wb'))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_seq_path",
        type=str,
        help="input probability sequence path"
    )
    parser.add_argument(
        "--input_prob_path",
        type=str,
        help="input probability sequence path"
    )
    parser.add_argument(
        "--output_path",
        type=str,
        help="output storage path",
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    generate_reference_seq_matrix(args.input_seq_path, args.input_prob_path, args.output_path)


if __name__ == '__main__':
    main()
