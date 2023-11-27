import sys
import numpy as np
import pandas as pd
import h5py
import os
import re

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <junction_tsv> <all_coords>\n' % sys.argv[0])
    sys.exit(1)
infname = sys.argv[1]
coordfname = sys.argv[2]
outfname = re.sub(r'.tsv.gz$', '', infname) + '.projected.hdf5'

### load coordinates of all introns
print(f'loading coords from {coordfname}')
coords = pd.read_csv(coordfname, sep='\t')

### load data
data = pd.read_csv(infname, sep='\t')

OUT = h5py.File(outfname, 'w')

for chrm in data['chr'].unique():
    for strand in ['+', '-']:
        d = data[(data['chr'] == chrm) & (data['strand'] == strand)]
        c = coords[(coords['chr'] == chrm) & (coords['strand'] == strand)]

        ### check sorting
        c = c.sort_values(by=['junction_start', 'junction_end'])
        d = d.sort_values(by=['junction_start', 'junction_end'])

        cs = np.array(c['junction_start'].astype('str') + ':' + c['junction_end'].astype('str'), dtype='str')
        ds = np.array(d['junction_start'].astype('str') + ':' + d['junction_end'].astype('str'), dtype='str')
        idx = np.where(np.in1d(cs, ds))[0]
        assert idx.shape[0] == ds.shape[0] ## make sure we find all junctions
    
        counts = np.zeros((cs.shape[0],), dtype='int')
        counts[idx] = d['count'].values
        
        print('writing %s %s' % (chrm, strand))
        OUT.create_dataset(name=chrm+':'+strand, data=counts, compression='gzip')
OUT.close()
        
