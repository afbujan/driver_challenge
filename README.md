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


The module can also be used from another module by using the approprite import:

```
import concatenate_reads as cr
concatenated_sequence = cr.getCompleteSequence(reads, threshold=3)
```

That would return a string with the concatenated DNA sequence.

