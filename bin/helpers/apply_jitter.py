import pandas as pd

import numpy as np

import sys



masterlist = pd.read_table(
    sys.argv[1],
    names=['chrom', 'start', 'end', 'dhs_id', 'summit']
)

sample_id = sys.argv[2]

jitters = np.load(sys.argv[3], mmap_mode='r')[sample_id]
dhs_idx = jitters[:, 0]
jitter_offsets = jitters[:, 1]

# keep only the DHSs sampled for this sample_id, drop the rest
sampled_masterlist = masterlist.iloc[dhs_idx].copy()
sampled_masterlist['summit'] += jitter_offsets.astype(int)

sampled_masterlist['start'] = sampled_masterlist['summit'].astype(int)
sampled_masterlist['end'] = sampled_masterlist['start'] + 1

sampled_masterlist = sampled_masterlist.sort_values(['chrom', 'start', 'end'])

sampled_masterlist[['chrom', 'start', 'end', 'dhs_id', 'summit']].to_csv(
    sys.argv[4],
    sep='\t',
    header=False,
    index=False,
)

# original masterlist row index for each row of the sorted bed above, for scattering results back
np.savetxt(sys.argv[5], sampled_masterlist.index.to_numpy(), fmt="%d")