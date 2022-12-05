build:
	gcc input_generator.c -o inmaker.out
	mpicc -g main.c

run:
	./inmaker.out 500 10000
	mpirun -n 2 a.out -i < data.in