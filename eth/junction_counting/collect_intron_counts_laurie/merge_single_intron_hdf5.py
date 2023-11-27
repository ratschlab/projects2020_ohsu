import sys
import numpy as np
import pandas as pd
import h5py
import glob
import os
import copy

if len(sys.argv) < 4:
    sys.stderr.write('Usage: %s <hdf5_dir> <outname> <chrm_list (comma-sep)>\n' % sys.argv[0])
    sys.exit(1)
indir = sys.argv[1]
outfname = sys.argv[2]
chrms = sys.argv[3].split(',')

### set up output structure
OUT = h5py.File(outfname, 'w')
OUT.create_dataset(name='chrms', data=np.ones(shape=(0,), dtype='|S10'), compression='gzip', maxshape=(None,))
OUT.create_dataset(name='pos', data=np.ones(shape=(0, 3), dtype='int'), compression='gzip', maxshape=(None, 3))

for chrm in chrms:
    print(f'Collecting {chrm}')
    fname = os.path.join(indir, f'all_junctions.projected.{chrm}.hdf5')
    IN = h5py.File(fname, 'r')
    if not 'counts' in OUT:
        OUT.create_dataset(name='counts', data=np.ones(shape=(0, IN['samples'].shape[0]), dtype='int'), compression='gzip', maxshape=(None, IN['samples'].shape[0]))
    if not 'samples' in OUT:
        OUT.create_dataset(name='samples', data=IN['samples'][:], compression='gzip')
    else:
        assert np.all(OUT['samples'][:] == IN['samples'][:])
    # collect
    print('...collect')
    print('......coords')
    starts = np.hstack([IN[f'{chrm}:+:junction_start'][:], IN[f'{chrm}:-:junction_start'][:]])
    ends = np.hstack([IN[f'{chrm}:+:junction_end'][:], IN[f'{chrm}:-:junction_end'][:]]) 
    strands = np.hstack([['+'] * IN[f'{chrm}:+:junction_start'].shape[0], ['-'] * IN[f'{chrm}:-:junction_start'].shape[0]]) 
    strands = [1 if _ == '+' else 2 for _ in strands]
    print('......counts')
    counts = np.vstack([IN[f'{chrm}:+:count'][:].astype(np.int16), IN[f'{chrm}:-:count'][:].astype(np.int16)])
    # aggregate
    print('...aggregate')
    coords = np.vstack([starts, ends, strands]).T
    del starts, ends, strands
    # sort
    s_idx = np.lexsort([coords[:, -i] for i in range(1, coords.shape[1] + 1)])
    coords = coords[s_idx, :]
    counts = counts[s_idx, :]
    # write
    print('...write')
    shp = OUT['pos'].shape
    OUT['pos'].resize((shp[0] + coords.shape[0], shp[1]))
    OUT['pos'][shp[0]:, :] = coords
    OUT['chrms'].resize((shp[0] + coords.shape[0],))
    OUT['chrms'][shp[0]:] = [chrm] * coords.shape[0]
    shp = OUT['counts'].shape
    OUT['counts'].resize((shp[0] + coords.shape[0], shp[1]))
    OUT['counts'][shp[0]:, :] = counts
    del counts
OUT.close()
