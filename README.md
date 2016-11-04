Driver Group challenge
=====================

Solution to the Driver coding challenge
---------------------------------------

The solution to the challenge is a module called ```concatenate_reads.py```

This python module can be run from the comman line by typing:

```
python concatenate_reads.py --help
``` 

This command will display the module's help info.


You can concatenate reads from a specific file by typing:


```
python concatenate_reads.py -f path/to/your/file.fasta
``` 

The output of the program will print the concatenated DNA sequence.


The module can also be used from another module with the appropriate import:

```
import concatenate_reads as cr
reads = cr.loadReads(filename)
concatenated_sequence = cr.getCompleteSequence(reads, threshold=3)
```

That would return a string with the concatenated DNA sequence.


How does the program work
-------------------------

First, this program assumes that there are no sequencing errors and that there
are no reads nested within reads.

These are the steps needed to solve the problem: 

 - Loading the reads and storing them in a dictionary
 - Creating a matrix of overlap lengths, which tells us the length of the overlap between all pairs of reads in both directions
 - Finding the order of the reads in the complete sequence, which involves: 
    - Finding the left-most read (the one with no significant** overlap to the left) 
    - Recursively finding the reads to the right of the subsequent read
 - Stiching all the reads together

**significant is in the program controlled by a threshold parameter. In general, an overlap larger than 3 base pairs could be considered significant since it is very unlikely to happen by chance. However, the threshold can be set higher if wanted. In any case, the way the program is implemeted the threshold is fixed for all reads.

Testing the program
-------------------

The program can be tested using the test module ```concatenate_reads_test.py```

