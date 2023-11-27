import h5py
import sys
import re
import numpy as np

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <spladder_junctions.hdf5>\n' % sys.argv[0])
    sys.exit(1)
infile = sys.argv[1]

outfile = re.sub(r'.hdf5$', '', infile) + '.tsv.gz'
data = [['chr', 'strand', 'junction_start', 'junction_end', 'count']]
with h5py.File(infile, 'r') as IN:
    for k in sorted(IN.keys()):
        if not (k.endswith('introns_p') or k.endswith('introns_m')):
            continue
        strand = '+' if k.endswith('introns_p') else '-'
        chrm = k.split('_')[0]
        sys.stderr.write('parsing %s\n' % k)
        for rec in IN[k][:]:
            data.append([chrm, strand, str(rec[0]), str(rec[1]), str(rec[2])])
sys.stderr.write('saving to %s\n' % outfile)
np.savetxt(outfile, np.array(data, dtype='str'), delimiter='\t', fmt='%s')


