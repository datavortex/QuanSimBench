### Please test your computer as follows:

1- First run the benchmark until you get the highest performance in the column "States/s". That usually corresponds to the 
largest number of nodes, or to the largest number of qubits that you can simulate (column Qubits). 
Verify that you get the "Pass=yes" result on all rows.

2- Submit the following fields in a text file to santiago.betelu@datavortex.com 
- Your name & institution
- Computer name
- CPU name
- Network type
- Memory per node
- Number of nodes used
- Physical cores per node used
- Total physical cores used
- Qubits (as said before, for which states/s is largest) 
- Rawtime 
- States/s 
- States/s/core
- Power usage (if known)


Example:
```
- Santiago Betelu, UNT
- UT's Stampede2
- Xeon Phi 7250
- Infiniband network
- 96GB/node
- 256 nodes used
- 64 cores used
- total cores 16384
- Qubits: 41
- States/s: 4.4e11
- States/s/core: 2.7e7
- Power usage unknown
```
