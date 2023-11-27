import sys
import numpy as np
import pandas as pd
import h5py
import glob
import os
import copy

if len(sys.argv) < 5:
    sys.stderr.write('Usage: %s <hdf5_dir> <chrm> <all_coords.tsv> <outname>\n' % sys.argv[0])
    sys.exit(1)
indir = sys.argv[1]
coordfname = sys.argv[3]
outfname = sys.argv[4]

### load coordinates of all introns
coords = pd.read_csv(coordfname, sep='\t', dtype='str')

if sys.argv[2] == 'all':
    chrms = coords['chr'].unique()
else:
    chrms = [sys.argv[2]]

### set up output structure
OUT = h5py.File(outfname, 'w')
strands = coords['strand'].unique()
for chrm in chrms:
    for strand in strands:
        c = coords[(coords['chr'] == chrm) & (coords['strand'] == strand)]
        c = c.sort_values(by=['junction_start', 'junction_end'])
        OUT.create_dataset(name=chrm+':'+strand+':count', data=np.ones(shape=(c.shape[0], 0), dtype='int'), compression='gzip', maxshape=(c.shape[0], None))
        OUT.create_dataset(name=chrm+':'+strand+':junction_start', data=c['junction_start'], dtype='int', compression='gzip')
        OUT.create_dataset(name=chrm+':'+strand+':junction_end', data=c['junction_end'], dtype='int', compression='gzip')

infiles = glob.glob(os.path.join(indir, '*.projected.hdf5'))
samples = []
for i, fname in enumerate(infiles):
    print(fname + ' (%i/%i)' % (i, len(infiles)))
    IN = h5py.File(fname, 'r')
    for chrm in chrms:
        for strand in strands:
            shp = OUT[chrm+':'+strand+':count'].shape
            OUT[chrm+':'+strand+':count'].resize((shp[0], shp[1] + 1))
            if chrm+':'+strand in IN:
                OUT[chrm+':'+strand+':count'][:, shp[1]] = IN[chrm+':'+strand][:]
            else:
                OUT[chrm+':'+strand+':count'][:, shp[1]] = np.zeros((shp[0], ), dtype='int')
    samples.append(os.path.basename(fname).split('.')[0])
OUT.create_dataset(name='samples', data=[_.encode('utf-8') for _ in samples])
OUT.close()
