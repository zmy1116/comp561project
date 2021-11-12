import numpy as np
import pickle
# %% Example 1

query1 = "AAAAAAAAAAAA"

ref1 = np.zeros((16, 4))
for k in range(4): ref1[k][0] = 1
for k in range(3, 5): ref1[k][3] = 1
for k in range(5, 13): ref1[k][0] = 1
for k in range(13, 15): ref1[k][3] = 1

pickle.dump(ref1, open("ref1.p", 'wb'))


"""
RUN
blast("ref1.p", query1, gap_penalty=-1, gap_bias=-2,
          hashtable_method='consensus_seed_seq',
          k=3, delta=1, mismatch_score=5, score_method='sum_proba_score',
          substitution=dict(), output_path="")
"""


# Expected result
# AAA--AAAAAAAA
# AAATTAAAAAAAA(TT)
#


# %% Example 2
query2 = "AAAAAAAAAAAAAA"

ref2 = np.zeros((15, 4))
for k in range(3): ref2[k][0] = 1
for k in range(3, 5): ref2[k][3] = 1
for k in range(5, 13): ref2[k][0] = 1
for k in range(13, 15): ref2[k][3] = 1

pickle.dump(ref2, open("ref2.p", 'wb'))

mismatch_score = -5

"""
RUN
blast("ref2.p", query2, gap_penalty=-1, gap_bias=-2,
          hashtable_method='consensus_seed_seq',
          k=3, delta=1, mismatch_score=5, score_method='sum_proba_score',
          substitution=dict(), output_path="")
"""

# Expected result
# AAA--AAAAAAAAAAA
# AAATTAAAAAAAA---TT
#
