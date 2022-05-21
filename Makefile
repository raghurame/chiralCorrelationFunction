all:
	gcc -o main main.c -lm
	./main chain.lammpstrj dihedral600.dump
