# 2d-fft

## Description of files

1. serial.c - Sequential version
2. parallel_simple.c - Parallel version using MPI_Send and MPI_Recv
3. parallel_collective.c - Parallel version using collective methods.
4. parallel_tasks.c - Parallel version using 2 processors for each task, has to be run with 8 processor.

## Compile

```bash
$ gcc -o serial serial.c -lm

$ mpicc -o parallel_simple ./parallel_simple.c -lm

$ mpicc -o parallel_collective ./parallel_collective.c -lm

$ mpicc -o parallel_tasks ./parallel_tasks.c -lm
```

## Running

```bash
$ ./serial

$ mpirun -n <Number of processors> ./parallel_simple

$ mpirun -n <Number of processors> ./parallel_collective

$ mpirun -n <Number of processors> ./parallel_tasks
```
