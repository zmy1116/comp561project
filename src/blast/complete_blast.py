import numpy as np
import pickle
import json
from src.utils.generate_query import generate_queries
from src.blast.hashtable_generation_inference import consensus_seq_match
from src.blast.ungapped_extension import ungapped_extension
from src.blast.gapped_extension_j_mingmodified import gapped_extension


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
