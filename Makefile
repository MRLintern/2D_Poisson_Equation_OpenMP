CC = gcc
CFLAGS = -g -O2 -Wall -Wextra -fopenmp


all: poisson_main
	
main: main.c
	$(CC) $(CFLAGS) -o poisson_main poisson_main.c -lm
	
clean:
	poisson_main
