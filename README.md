# QuanSimBench

## Why Quantum factorization simulation as a benchmark for HPC
- Quantum factorization a widely known & relevant problem
- The goal of this test is to quantify the ability of a computer to simulate ideal quantum circuits
- Easy to validate and understand the results
- Each additional qubit doubles RAM usage, CPU power and internode communication: good to test large machines
- Runs in a reasonable time: 1-3 hours
- Portable with less than 300 lines of C and MPI
- It just runs: no input or special knowledge from user
- Runs from a laptop to a large supercomputer

## How to compile
```
mpicc -Ofast quansimbench.c -o quansimbench -lm -Wall
sbatch quansimbench.batch
```
You may want to add optimization flags for your architecture

## Batch file quansimbench.batch
```
#SBATCH –o outputfile
#SBATCH --nodes=8     # 8 nodes
#SBATCH –n 256        # 64 cores per node
#SBATCH –p normal     # for KNL processors
#SBATCH –t 02:00:00   # usually less than 3 hours
ibrun ./quansimbench
```
