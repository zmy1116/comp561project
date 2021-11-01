# comp561project
project for comp561


```buildoutcfg
export PYTHONPAHT= LOCAL_PACKAGE_PATH 
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
python /home/mcb/users/mzhou3/comp561project/src/preprocessing.py \
	--input_seq_path=/home/mcb/users/mzhou3/comp561project/src/data/chr22.maf.ancestors.42000000.complete.boreo.fa \
	--input_prob_path=/home/mcb/users/mzhou3/comp561project/src/data/chr22.maf.ancestors.42000000.complete.boreo.conf \
	--output_path=/home/mcb/users/mzhou3/comp561project/src/data/reference_matrix.p
```

