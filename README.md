# DNA Sequence Alignment

The project provides two solutions for the DNA sequence alignment problem:

- Basic Version using Dynamic Programming having the time complexity of O(mn) and space complexity of O(mn).
- Efficient Version using Dynamic Programming and Divide and Conquer having the time complexity of O(mn) and space complexity of O(n). 

## How to run
- Install the dependencies in requirement.txt
- Generate input file for the program as shown in the [doc](https://github.com/kishanmurthy/dna-sequence-alignment/blob/main/docs/CSCI570_Fall2021_FinalProject.pdf).
- Run the basic version (basic.py) and efficient version (efficient.py) with input file as argument.
- Example
  1. ```python basic.py data/input1.txt```
  2. ```python efficient.py data/input1.txt```

## Summary
After running both the basic and efficient solutions for the DNA Sequence Alignment for the problem size from 2 to 200 and have made the following observations:

<img src="https://github.com/kishanmurthy/dna-sequence-alignment/blob/main/outputs/CPUPlot.png" width="700">

**1. Comparing CPU time with increasing Problem Size**
   - Initially, for problem size of length 2, the time taken by the basic version of finding an optimal DNA sequence alignment is 0.02 milliseconds, and the efficient version of finding an optimal DNA sequence alignment takes 0.04 milliseconds.
   - When the problem size increases to 100, we observe that the time taken by an efficient version of the algorithm is 0.021 seconds, while the time taken by the basic version of the algorithm is 0.008 seconds.
   - When the problem size increases to 200, we observe that the time taken by an efficient version of the algorithm is 0.085 seconds, while the time taken by the basic version of the algorithm is 0.034 seconds. 
   - The general trend observed is that the time taken for the solution is at a constant factor between basic and efficient algorithms. We can generalize and say that time taken by an efficient algorithm takes 2.5x the time taken by the basic algorithm given the problem size.

<img src="https://github.com/kishanmurthy/dna-sequence-alignment/blob/main/outputs/MemoryPlot.png" width="700">

**2. Comparing Memory Usage with increasing Problem Size**
   - Initially, for problem size of length 2, the memory usage for the basic version of finding the sequence alignment is 0.672 KB, which is lower than the efficient version of finding the sequence alignment at 0.768 KB.
   - When the problem size increases to 100, the memory usage explodes for the basic version of the algorithm, taking 90.79 KB, whereas the efficient version takes 6.17 KB. There is a 14x times difference between the basic and memory-efficient versions.
   - When the problem size increases to 200, the memory usage explodes further, and there is a 30x difference between memory consumption of basic and efficient algorithms. The basic version takes 372.61 KB, while the memory-efficient version takes 12.07 KB only.
   - The general trend observed is that total memory required to return an optimal sequence alignment increases squared polynomially with increasing problem size for the basic version of the code. The memory consumption of the basic version is infeasible for the human genome sequence matching, given the problem size in billions.
   - The efficient version of finding an optimal DNA sequence alignment performs better than the basic version, and it increases linearly with the increasing problem size.


**The combination of dynamic programming and divide and conquer provides a scalable and efficient solution to the DNA Sequence alignment problem. This makes the efficient solution apt for human genome sequence matching.**
