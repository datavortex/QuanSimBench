# QuanSimBench

## Why Quantum factorization simulation as a benchmark for HPC
- Quantum factorization a widely known & relevant problem
- The goal is to quantify the ability of a computer to simulate ideal quantum circuits
- Easy to validate and understand the results
- Each additional qubit doubles RAM usage, CPU power and internode communication: good to test large machines
- Runs in a reasonable time: 30'-3 hours
- Portable with less than 300 lines of C and MPI
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
#SBATCH –o outputfile-%J
#SBATCH --nodes=8     # 8 nodes (must be a power of 2)
#SBATCH –n 256        # 64 cores per node (must be a power of 2)
#SBATCH –p normal     # what queue
#SBATCH –t 02:00:00   # usually less than 3 hours
ibrun ./quansimbench
```
