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
```
gcc -Ofast quansimbench-sharemem.c -o quansimbench -lm -Wall -fopenmp``
export OMP_NUM_THREADS=32
./quansimbench
```
