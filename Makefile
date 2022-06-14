CC = gcc
CFLAGS = -g -O2 -Wall -Wextra -fopenmp


all: poisson
	
poisson: poisson.c
	$(CC) $(CFLAGS) -o poisson poisson.c -lm
	
clean:
	poisson
