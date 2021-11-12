import src.blast.hashtable_generation_inference as hash
import src.blast.ungapped_extension as ungapped
import src.blast.gapped_extension as gapped
# %%
def blast(reference_matrix_file, query, gap_penalty=-1, gap_bias=-2,
          hashtable_method='consensus_seed_seq',
          k=11, delta=0.1, score_method='sum_proba_score', substitution=dict(),
          output_path=""):

    table_data = hash.hashtable_generation(reference_matrix_file,
                                           hashtable_method, k,
                                           output_path)

    print(1)
    matches_dict = hash.consensus_seq_match(table_data, query)

    print(2)
    ungapped_extensions = ungapped.ungapped_extension(query, matches_dict,
                                                      reference_matrix_file, k,
                                                      delta, score_method,
                                                      substitution)

    print(3)

    gapped_extensions = gapped.gapped_extension(query, reference_matrix_file,
                                                ungapped_extensions,
                                                score_method, substitution,
                                                gap_penalty, gap_bias)
    print(4)
    return gapped_extensions


# %%
import os
os.chdir("Documents/McGill_courses/COMP561/group_project/comp561project")

# %%
reference_matrix_file = "src/data/reference_matrix.p"
query = "A" * 26
blast(reference_matrix_file, query, gap_penalty=-1, gap_bias=-2,
          hashtable_method='consensus_seed_seq',
          k=11, delta=0.1, score_method='sum_proba_score', substitution=dict(),
          output_path="")