import numpy as np
from Bio import SeqIO
import pdb
from optparse import OptionParser

def main():
    usage = '%prog [options]'
    parser = OptionParser(usage)
    parser.add_option("-f",'--file_path',type="string",default='',\
        help="path to the DNA file to be concatenated.")
    (options, args) = parser.parse_args()
    printConcatReads(options.file_path)

def printConcatReads(file_path):
    """Returns the concatenate sequence of DNA reads
    Args:
       file_path: path to the DNA file
    """
    reads = loadReads(file_path)
    overlapLengths = getOverlapLengths(reads)
    order = getOrder(getLeftMostID(overlapLengths), overlapLengths)
    print '\nFor file: %s'%file_path
    print '\nThe complete DNA sequence is: [%s]'%concatReads(order, reads, overlapLengths)

def get_records(file_path, file_format='fasta'):
    """Loads a DNA file
    Args:
        file_path (string): path to the DNA file
        file_format (string): format of the DNA file
    Returns:
        records (list): DNA reads
    """
    records = SeqIO.parse(open(file_path, 'rU'), file_format)
    records = list(records)
    return records

def loadReads(filename, fileformat='fasta'):
    """Loads DNA data from file
    Args:
        filename (string): name of the file
        fileformat (string): format of the file
    Returns:
        reads (dict): a dictionary with the DNA reads
    """
    records = get_records(filename, fileformat)
    records = [list(r) for r in records]
    reads = {}
    for i, r in enumerate(records):
        reads[str(i)]=''.join(r)
    return reads

def getOverlap(seq_left, seq_right):
    """Returns the overlapping section between 2 DNA seqs
    Args:
        seq_left (string): left sequence
        seq_right (string): right sequence
    Returns:
        overlap (string): overlapping sequence
    """
    overlap = ''
    seq_len = len(seq_left)
    for i in xrange(seq_len):
        if seq_left[i:] == seq_right[:seq_len-i]:
            overlap = seq_left[i:]
            break
    return overlap

def getOverlapLengths(reads):
    """Computes the length of all the overlaps between pairs of DNA seqs
        from left-right (rows) and right-left (columns)
    Args:
        reads (dict): all the DNA reads to be compared
    Returns:
       overlapLengths (ndarray): the length of all pairwise overlaps 
    """
    n_reads = len(reads)
    overlapLengths = np.zeros((n_reads, n_reads), dtype='int')
    for i, s1 in enumerate(reads.values()):
        for j, s2 in enumerate(reads.values()):
            if i==j:
                continue
            overlapLengths[i,j] = len(getOverlap(s1, s2))
    return overlapLengths

def getLeftMostID(overlapLengths, threshold=300):
    """Returns the id of the starting read of the sequence, ie: the 
        first read that has no significant (>threshold) overlapping sequences 
        to the left (columns) but has a significant overlap to the right.
    Args:
        overlapLegths (ndarray): length of pairwise overlaps between reads
        threshold (int): number of bases that makes an overlap significant
    Returns:
        leftMost_id (int): index of the left-most read
    """
    leftMost_id = None
    for i in xrange(overlapLengths.shape[0]):
        if not np.any(overlapLengths[:,i]>threshold) and np.any(overlapLengths[i]>threshold):
            leftMost_id = i
            break
    assert leftMost_id is not None, 'Left-most sequence could not be found!'
    return leftMost_id

def getOrder(idx, overlapLengths, threshold=300):
    """Returns the sequence order
    Args:
        idx (int): the id of the left most read in the sequence
        overlapLegths (ndarray): length of pairwise overlaps between reads
        threshold (int): number of bases that makes an overlap significant
    Returns:    
        order (list): ordered indexes of DNA fragments
    """
    if overlapLengths[idx].max() < threshold:
        return [idx]
    else:
        nextRead = np.argmax(overlapLengths[idx])
        return [idx] + getOrder(nextRead, overlapLengths)

def concatReads(order, reads, overlapLengths):
    """Returns a concatenated sequence of DNA reads
    Args:
        order (list): ordered indexes of DNA reads
        reads (dict): all the DNA reads
        overlapLegths (ndarray): length of pairwise overlaps between DNA reads
    Returns:    
        completeSequence (string): complete DNA sequence
    """ 
    completeSequence = ''
    for i in xrange(len(order)-1):
        rightOverlap = overlapLengths[order[i],order[i+1]]
        completeSequence += reads[str(order[i])][:-rightOverlap]
    completeSequence += reads[str(order[-1])]
    return completeSequence


if __name__ == "__main__":
    main()
