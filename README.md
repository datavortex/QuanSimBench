## What this benchmark does

- A gate-by-gate simulation of a quantum computer with a full state vector formulation. 
- It simulates a simplified version of Shor's algorithm with increasing number of qubits until resources are exhausted.
- The score (states/s) is how fast the Approximate Quantum Fourier Transform can be computed (the modular exponentiation is not timed).  
- The goal is to quantify the ability of a computer to simulate ideal quantum circuits.

## Why quantum factorization simulation as a benchmark for HPC
- Quantum factorization is a widely known & relevant problem, and it is known to be hard to simulate.
- Each additional qubit doubles RAM usage, CPU power and internode communication: good to test large machines.
- The output is reproducible, easy to validate and understand.
- Runs in a reasonable time: 1/2‒3 hours.
- Minimal portable code with less than ~300 lines of C and MPI, very few dependencies.
- It just runs: no input or specialized knowledge from user.

## How to compile and run
This is system dependent. For example, in systems with mpicc and Slurm one usually compiles with
```bash
mpicc -Ofast quansimbench.c -o quansimbench -lm
```
and runs the program with
```bash
sbatch quansimbench.batch
```
using a Slurm batch script such as:
```bash
#SBATCH –o output-%J.txt
#SBATCH --nodes=8     # 8 nodes (must be a power of 2)
#SBATCH –n 64         # 8 ranks per node (must be a power of 2)
#SBATCH –p normal     # what queue
#SBATCH –t 02:00:00   # usually less than 3 hours
mpirun ./quansimbench
```
**In all systems, the number of nodes and number of cores must be a power of two.**

The code contains comments with hints to optional addition of OpenMP as well.

## How to test the ability of your computer to simulate quantum AQFTs

1. First run the benchmark until completion with increasing number of nodes/ranks until you get the largest number of States/s in the last row. That usually corresponds to the largest number of nodes/ranks you can run. If the program hangs, report the last visible row.

2. Verify that you get the "Pass=yes" result on all rows.

3. Submit your results to [QuanSimBench Submission](https://docs.google.com/forms/d/e/1FAIpQLSeVwp_4FZJWyS5UsfBrtxq8PXkKJLoRvgHkpfTuOuJ-wcudiw/viewform?usp=sf_link).

4. The results will soon be posted on our GitHub site.
