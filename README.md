Driver Group challenge
=====================

Solution to the Driver Group coding challenge
---------------------------------------

The solution to the challenge is a module called ```concatenate_reads.py```. This python module can be run from the command line by typing:

```
python concatenate_reads.py --help
``` 

The command above will display the module's help info. You can concatenate reads from a specific file by typing:

```
python concatenate_reads.py -f path/to/your/file.fasta
``` 

The output of the program will print the concatenated DNA sequence. The module can also be used within another module with the appropriate import statement:

```
import concatenate_reads as cr
reads = cr.loadReads(filename)
concatenated_sequence = cr.getCompleteSequence(reads, threshold=3)
```

These lines will return a string with the concatenated DNA sequence.


How the program works
---------------------

The program assumes that there are no sequencing errors and that there are no reads nested within reads. The steps needed to solve the problem are the following: 

 - Loading the reads and storing them in a dictionary
 - Creating a matrix of overlap lengths, which contains the length of the overlap between all pairs of reads in both directions
 - Finding the order of the reads in the complete sequence, which involves: 
    - Finding the left-most read (the one with no significant** overlap to the left) 
    - Recursively finding the read with significant overlap to the right of the subsequent read
 - Stiching all the reads together once the order is known

**significant is (in the program) controlled by a threshold parameter. In general, an overlap larger than 3 base pairs could be considered significant since it is very unlikely to happen by chance. However, the threshold can be set higher if needed. In the program, the threshold is always the same across reads.

Testing the program
-------------------

The program can be tested using the test module ```concatenate_reads_test.py```. This module can be run from the command line by typing:

```
python concatenate_reads_test.py --help
```

The command above will display the module's help info. There are 3 parameters that can be controlled: the number of reads, the length of the reads, and the ratio between the length of the overlapping and non-overlapping regions (called beta). To test the code you can type the following:

```
python concatenate_reads_test.py --n 50 --l 1000 
```

If no parameters are given the default ones are used. The test module generates a series of overlapping reads and mixes them randomly. The output of ```concatenate_reads.py``` is then compared with the ground truth. The output of ```concatenate_reads_test.py``` will print to the screen whether the test passed, a list with the paramenters of the test, and the time it took for the algorithm to complete the task. An example output is shown below:

```
Testing with following parameters:
	*n_reads:     50
	*read_length: 1000
	*beta:        0.2
Test passed:  True  elapsed time: 1.1814 seconds
```
