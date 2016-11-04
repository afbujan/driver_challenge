import numpy as np
import concatenate_reads as cr
import pdb, time
from optparse import OptionParser

def main():
    usage = '%prog [options]'
    parser = OptionParser(usage)
    parser.add_option("-l",'--read_length',type="int",default=1000,\
        help="lenght of the DNA reads")
    parser.add_option("-n","--n_reads",type="int",default=50,\
        help="number of DNA reads")
    parser.add_option("-b","--beta",type="float",default=.5,\
        help="length ratio between overlapping and non-overlapping regions")
    (options, args) = parser.parse_args()
    beta = options.beta
    assert beta>=0. and beta<=1., 'Beta ratio needs to be a number between 0 and 1!'
    test(n_reads=options.n_reads,read_length=options.read_length,beta=beta)

def test(n_reads=50, read_length=1000, beta=.5):
    reads = {}
    idx = np.random.permutation(np.arange(n_reads))
    overlapping_length = int(np.floor(read_length*beta))
    non_overlapping_length = read_length-overlapping_length
    overlapping = np.random.randint(10,size=(n_reads,overlapping_length)).astype('S8')
    groundTruth = ''
    for i in xrange(n_reads):
	non_overlapping = ''.join((non_overlapping_length)*[str(i)])
	if i==0:
	    reads[str(idx[i])] = ''.join(non_overlapping)+''.join(overlapping[i].tolist())
	    groundTruth += ''.join(non_overlapping)
	elif i==n_reads-1:
	    reads[str(idx[i])] = ''.join(overlapping[i-1].tolist())+''.join(non_overlapping)
	    groundTruth += reads[str(idx[i])]
	else:
	    reads[str(idx[i])] = ''.join(overlapping[i-1].tolist())+''.join(non_overlapping)+''.join(overlapping[i].tolist())
	    groundTruth += ''.join(overlapping[i-1].tolist())+''.join(non_overlapping) 
    tic = time.time()
    print "Testing with following parameters:"
    print "\t*n_reads:     %i"%n_reads
    print "\t*read_length: %i"%read_length
    print "\t*beta:        %.1f"%beta
    completeSequence = cr.getCompleteSequence(reads, threshold=overlapping_length-1)
    print 'Test passed: ', completeSequence==groundTruth,' elapsed time: %.4f seconds'%(time.time()-tic)

if __name__ == "__main__":
    main()
