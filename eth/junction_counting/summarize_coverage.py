import sys

chrm = None
cnt = 0
cov = 0
for line in sys.stdin:
    sl = line.strip().split('\t')
    if sl[0] != chrm:
        if not chrm is None:
            print('\t'.join([chrm, str(cnt), str(cov)]))
        chrm = sl[0]
        cnt = 0
        cov = 0
    cnt += 1
    cov += int(sl[1])
if not chrm is None and cnt > 0:
    print('\t'.join([chrm, str(cnt), str(cov)]))
