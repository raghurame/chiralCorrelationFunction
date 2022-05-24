all:
	export OMP_DYNAMIC=true
	gcc -c chiralCorrelation.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c main.c -Wall -fstack-protector -g -fopenmp -lm
	gcc main.o chiralCorrelation.o -o main -lm -fopenmp -Wall -g -fstack-protector
	./main chain.lammpstrj dihedral600.dump
