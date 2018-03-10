gcc ga-ac-seq.c -o seq -lm

mpicc -fopenmp ga-ac-parallel.c -o par -lm

