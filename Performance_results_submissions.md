### Please test your computer as follows:

A) First run the benchmark with increasing number of nodes/cores until you get the highest performance in the column "States/s". That usually corresponds to the largest number of qubits that you can simulate (column Qubits), and the largest number of nodes you can run. Consider that a simulation with Q qubits will need 8x2^Q bytes of global memory. 

B) Verify that you get the "Pass=yes" result on all rows.

C) Select the row with the highest performance and submit the following fields in a text file to santiago.betelu@datavortex.com 
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
- Attach a text file "output.txt" with the output of the program.

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
