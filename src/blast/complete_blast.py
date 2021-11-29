import os
import pickle
import json
import argparse
from src.blast.hashtable_generation_inference import consensus_seq_match
from src.blast.ungapped_extension import ungapped_extension
from src.blast.gapped_extension_j_mingmodified import gapped_extension

DEFAULT_BLAST_ARGS = json.load(open(os.path.join(os.path.dirname(__file__), 'config', 'default_blast_args.json'), 'r'))


# %%
def blast(query, reference_matrix_file, seed_matching_args, ungapped_extension_args, gapped_extension_args):
    # step 1 seed matching

    # load seed matching table data
    table_path = seed_matching_args['seed_table_file']
    table_data = pickle.load(open(table_path, 'rb'))

    outputs_step1 = consensus_seq_match(table_data, query)

    # step 2 ungapped matching
    k = seed_matching_args['k']
    delta = ungapped_extension_args['delta']
    score_method = ungapped_extension_args['nt_score_method']
    mismatch_score = ungapped_extension_args['mismatch_score']
    substitution = ungapped_extension_args['substitution']
    outputs_step2 = ungapped_extension(query, outputs_step1, reference_matrix_file, k, delta,
                                       score_method, mismatch_score=mismatch_score, substitution=substitution)

    # step 3 gapped matching
    score_method = gapped_extension_args['nt_score_method']
    substitution = gapped_extension_args['substitution']
    gap_penalty = gapped_extension_args['gap_penalty']
    gap_bias = gapped_extension_args['gap_bias']
    mismatch_score = gapped_extension_args['mismatch_score']
    ref_max_length_factor = gapped_extension_args['ref_max_length_factor']
    outputs_step3 = gapped_extension(query, reference_matrix_file, outputs_step2, score_method,
                                     substitution, gap_penalty, gap_bias, mismatch_score, ref_max_length_factor)

    outputs_step3 = sorted(outputs_step3, key=lambda x: x['final_result']['score'] * -1)

    return outputs_step3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_file",
        type=str,
        help="output storage path"
    )

    parser.add_argument(
        "--reference_matrix_file",
        type=str,
        help="file path of the reference matrix file"
    )

    parser.add_argument(
        "--query_file",
        type=str,
        help="file path of the query, if None use a randomly generate query",
        default=None
    )

    parser.add_argument(
        "--config_file",
        type=str,
        help="file path of the algorithm configs, if None, use the default blast arguments file",
        default=None)
    parser.add_argument(
        "--hashtable_data_file",
        type=str,
        help="file path of the hashtable file, if None, use the default path",
        default=None
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.query_file is not None:
        query = open(args.query_file).readlines()[0]
    else:
        # use an example query
        query = 'TAAACCCCGGGCCATATTATGCTGGGAGGCTCACCTCTGTAATGTCGGGGTCTTGGAGGCTGAGGAGGTATTGCTTTGGGGAGTTCAAGCCCAGCCTGGG'

    if args.config_file is not None:
        blast_configs = json.load(open(args.config_file, 'r'))
    else:
        # use default args
        blast_configs = DEFAULT_BLAST_ARGS

    if args.hashtable_data_file is not None:
        args['seed_matching_args']['seed_table_file'] = args.hashtable_data_file

    reference_matrix_file = args.reference_matrix_file

    outputs = blast(query, reference_matrix_file, **blast_configs)
    for idx, out in enumerate(outputs):
        print('{0}  ref_start_idx:{1},   ref_end_idx:{2},    score: {3}'.format(idx,
                                                                                out['final_result']['ref_left_idx'],
                                                                                out['final_result']['ref_right_idx'],
                                                                                out['final_result']['score']))
    result = {
        'query': query,
        'result': outputs
    }
    output_file = args.output_file
    pickle.dump(result, open(output_file, 'wb'))


if __name__ == '__main__':
    main()
