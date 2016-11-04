import numpy as np
import concatenate_reads as cr
import pdb

n_reads = 50
read_length = 1000
reads = {}
idx = np.random.permutation(np.arange(n_reads))
overlapping = np.random.randint(10,size=(n_reads,read_length//2)).astype('S8')
for i in xrange(n_reads):
    non_overlapping = ''.join((read_length//4)*[str(i)])
    if i==0:
        reads[str(idx[i])] = ''.join(non_overlapping)+''.join(overlapping[i].tolist())
    elif i==n_reads-1:
        reads[str(idx[i])] = ''.join(overlapping[i-1].tolist())+''.join(non_overlapping)
    else:
        reads[str(idx[i])] = ''.join(overlapping[i-1].tolist())+''.join(non_overlapping)+''.join(overlapping[i].tolist())
cr.printConcatReads(reads)


