### Please test your computer as follows:

A) First run the benchmark with increasing number of nodes/cores until you get the highest performance in the column "States/s". That usually corresponds to the largest number of qubits that you can simulate (column Qubits), and the largest number of nodes you can run.  

B) Verify that you get the "Pass=yes" result on all rows.

C) Select the row with the highest performance and submit your results to <a href="https://docs.google.com/forms/d/e/1FAIpQLSeVwp_4FZJWyS5UsfBrtxq8PXkKJLoRvgHkpfTuOuJ-wcudiw/viewform?usp=sf_link" target="_blank">QuanSimBench Submission</a>


the following fields in a text file to santiago.betelu@datavortex.com 
- Your name & institution
- Computer name and vendor
- CPU name
- Network type
- Memory per node
- Year
- Number of nodes used (from batch file)
- Physical cores per node used 
- Total physical cores used (from batch file)
- Qubits (as said before, for which states/s is largest) 
- States/s 
- States/s/core
- Power usage (if known)
- Attach a text file "output.txt" with the output of the program.

Example:
```
- Santiago Betelu, UNT
- UT's Stampede2
- Xeon Phi 7250
- Fdr Infiniband network
- 96GB/node
- 2018
- 256 nodes used
- 64 cores used
- total cores 16384
- Qubits: 41
- States/s: 4.4e11
- States/s/core: 2.7e7
- Power usage unknown
```
