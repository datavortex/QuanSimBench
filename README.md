# QuanSimBench

- The test performs a gate-by-gate simulation of a quantum computer with a full state vector description. 
- It simulates a simplified version of Shor's algorithm with an increasing number of qubits Q until resources are exhausted 
- Only the Approximate Quantum Fourier Transform is timed, not the modular exponentiation

## Why Quantum factorization simulation as a benchmark for HPC
- Quantum factorization a widely known & relevant problem
- The goal is to quantify the ability of a computer to simulate ideal quantum circuits
- Easy to validate and understand the results
- Each additional qubit doubles RAM usage, CPU power and internode communication: good to test large machines
- Runs in a reasonable time: 1/2 - 3 hours
- Portable with less than ~300 lines of C and MPI, few dependencies
- It just runs: no input or special knowledge from user
- This reference implementation runs from a laptop to a large supercomputer

## How to compile
```
mpicc -Ofast quansimbench.c -o quansimbench -lm -Wall
sbatch quansimbench.batch
```
You may want to add optimization flags for your architecture

## Batch file quansimbench.batch
The number of nodes and number of cores must be a power of two:
```
#SBATCH –o output-%J.txt
#SBATCH --nodes=8     # 8 nodes (must be a power of 2)
#SBATCH –n 256        # 64 cores per node (must be a power of 2)
#SBATCH –p normal     # what queue
#SBATCH –t 02:00:00   # usually less than 3 hours
ibrun ./quansimbench
```
## Please test the ability of your computer to simulate quantum AQFT's as follows:

A) First run the benchmark with increasing number of nodes/cores until you get the highest performance in the column "States/s". That usually corresponds to the largest number of qubits that you can simulate (column Qubits), and the largest number of nodes you can run.  

B) Verify that you get the "Pass=yes" result on all rows.

C) Select the row with the highest performance and submit your results to <a href="https://docs.google.com/forms/d/e/1FAIpQLSeVwp_4FZJWyS5UsfBrtxq8PXkKJLoRvgHkpfTuOuJ-wcudiw/viewform?usp=sf_link" target="_blank">QuanSimBench Submission</a>

D) Soon you may see results posted in our github
