import src.blast.hashtable_generation_inference as hash
import src.blast.ungapped_extension as ungapped
import src.blast.gapped_extension as gapped

# %%
def blast(reference_matrix_file, query, hashtable_method, score_method, k=11,
          delta=0.1, score_method='sum_proba_score', substitution=dict()):

    table_data = hash.hashtable_generation(reference_matrix_file, method, k, output_path)
    matches_dict = consensus_seq_match(table_data, query)

    ungapped_extensions = ungapped.ungapped_extension(query, matches_dict,
                                                      reference_matrix_file, k,
                                                      delta, score_method,
                                                      substitution)

