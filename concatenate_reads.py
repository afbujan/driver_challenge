import numpy as np
from Bio import SeqIO
import pdb
from optparse import OptionParser

def main():
    usage = '%prog [options]'
    parser = OptionParser(usage)
    parser.add_option("-f",'--file_path',type="string",default='',\
        help="path to the DNA file to be concatenated.")
    parser.add_option("-t","--threshold",type="int",default=300,\
        help="min number of bases that make an overlap region significant")
    (options, args) = parser.parse_args()
    reads = loadReads(filename=options.file_path)
    print '\nDNA file: %s'%options.file_path
    printConcatReads(reads=reads,threshold=options.threshold)

def printConcatReads(reads,threshold=300):
    """Returns the concatenate sequence of DNA reads
    Args:
        threshold (int): number of bases that makes an overlap significant
        reads (dict): a dictionary with the DNA reads
    """
    overlapLengths = getOverlapLengths(reads)
    
    order = getOrder(getLeftMostID(overlapLengths, threshold=threshold), 
                     overlapLengths, threshold=threshold)
    print '\nComplete DNA sequence after concatenation: [%s]'%concatReads(order, reads, overlapLengths)

def getRecords(file_path, file_format='fasta'):
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
    """Loads reads of DNA from a file
    Args:
        filename (string): name of the file
        fileformat (string): format of the file
    Returns:
        reads (dict): a dictionary with the DNA reads
    """
    records = getRecords(filename, fileformat)
    records = [list(r) for r in records]
    reads = {}
    for i, r in enumerate(records):
        reads[str(i)]=''.join(r)
    return reads

def getReadLength(reads):
    """Returns the length of all the reads 
    Args:
        reads (dict): a dictionary with the DNA reads
    Returns:
        read_length (list): length of the DNA reads
    """
    read_length = []
    for i in reads:
        read_legth.append(len(reads[i]))
    return read_length

def getOverlap(left_read, right_read):
    """Returns the overlapping section between 2 DNA reads
    Args:
        left_read (string): left DNA read
        right_read (string): right DNA read
    Returns:
        overlap (string): overlapping DNA sequence
    """
    overlap = ''
    read_length = len(left_read)
    for i in xrange(read_length):
        if left_read[i:] == right_read[:read_length-i]:
            overlap = left_read[i:]
            break
    return overlap

def getOverlapLengths(reads):
    """Computes the length of all the overlaps between pairs of DNA reads 
        from left-right (rows) and right-left (columns)
    Args:
        reads (dict): all the DNA reads to be compared
    Returns:
       overlapLengths (ndarray): the length of all pairwise overlaps 
    """
    n_reads = len(reads)
    overlapLengths = np.zeros((n_reads, n_reads), dtype='int')
    for i, r1 in reads.items():
        for j, r2 in reads.items():
            if i==j:
                continue
            overlapLengths[int(i),int(j)] = len(getOverlap(r1, r2))
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
        overlapLength = overlapLengths[order[i],order[i+1]]
        completeSequence += reads[str(order[i])][:-overlapLength]
    completeSequence += reads[str(order[-1])]
    return completeSequence

if __name__ == "__main__":
    main()
