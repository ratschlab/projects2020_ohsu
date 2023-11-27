import sys
import h5py
import re

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <alignment.hdf5>\n' % sys.argv[0])
    sys.exit(1)
fname = sys.argv[1]
outfname = re.sub(r'.hdf5$', '', fname) + '.coverage_stat.tsv'

IN = h5py.File(fname, 'r')

keys = IN.keys()

cnt = 0
cov = 0
with open(outfname, 'w') as out:
    for k in sorted(keys):
        if not k.endswith('_reads_dat'):
            continue
        _cnt = IN[k].shape[0]
        _cov = int(IN[k][:].sum())
        out.write('\t'.join([k.split('_')[0], str(_cnt), str(_cov)]) + '\n')
        out.flush()
        cnt += _cnt
        cov += _cov
    out.write('\t'.join(['TOTAL', str(cnt), str(cov)]) + '\n')
    out.flush()
