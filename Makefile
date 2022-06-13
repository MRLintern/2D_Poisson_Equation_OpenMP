CC = gcc
CFLAGS = -g -O2 -Wall -Wextra -fopenmp


all: main
	
main: main.c
	$(CC) $(CFLAGS) -o main main.c -lm
	
clean:
	main
