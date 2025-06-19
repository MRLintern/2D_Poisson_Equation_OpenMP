# Makefile for Poisson Solver Project

# Compiler and flags
CC = gcc
CFLAGS = -O2 -fopenmp -Wall -std=c99

# Executable name
BIN = main

# Output directory for CSV files
RESULTS_DIR = results

# Source files
SRC = main.c PSolver.c

# Object files (optional)
OBJ = $(SRC:.c=.o)

# Default target: build the executable
all: $(BIN)

# Compile the program
$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $(BIN) $(SRC) -lm

# Run the program and ensure the results directory exists
run: $(BIN)
	mkdir -p $(RESULTS_DIR)     # Create the results directory if missing
	./$(BIN)

# Run unit tests using a Bash script
test:
	bash test.sh

# Clean up build artifacts and output files
clean:
	rm -f $(BIN)
	rm -rf $(RESULTS_DIR)
