# comp561project
project for comp561


```buildoutcfg
git clone https://github.com/zmy1116/comp561project
export PYTHONPATH = LOCAL_PACKAGE_PATH 
cd LOCAL_PACKAGE_PATH
```

## Data
- `src/data/chr22.maf.ancestors.42000000.complete.boreo.conf`: input probs data 
- `src/data/chr22.maf.ancestors.42000000.complete.boreo.fa`: input sequence data
-  `src/data/reference_matrix.p`: sequence probability matrix
```buildoutcfg
import pickle 
refernce_seq = pickle.load(open('src/data/reference_matrix.p', 'rb'))
```

## Modules 

- `preprocessing.py`  generate reference sequence matrix from raw data
```buildoutcfg
python src/preprocessing.py \
	--input_seq_path=src/data/chr22.maf.ancestors.42000000.complete.boreo.fa \
	--input_prob_path=src/data/chr22.maf.ancestors.42000000.complete.boreo.conf \
	--output_path=src/data/reference_matrix.p
```

- `complete_blast.py` do full blast 
```
python src/blast/complete_blast.py \
    --output_file=test \
    --reference_matrix_file=src/data/reference_matrix.p \
    --query_file=src/data/example_query.fa

```

## Notebooks
 
Two notebooks included in the repository, most of plots in our report are generated in these 2 notebooks:
-  `blast_evaluation_first_2_parts.ipynb`: test different hyperparameters/configuraiton for first 2 components of blast algorithm
-  `full_blast_evaluation.ipynb`: test the full blast algorithm in 3 set of queries:
    - 400 queries with no indels/substitutions of 4 difficulties, 100 length query
    - 400 queries with indels/substitutions of 4 difficulties, 100 length query 
    - 400 queries with indels/substitutions of 4 difficulties, 500 length query