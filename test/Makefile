# Compiler & compiler flags
CC = g++
CFLAGS = -Wall -lm

TEST_FP_FLAGS = \
    -fno-fast-math \
    -fno-unsafe-math-optimizations \
    -ffp-contract=off \
    -frounding-math \
    -fno-associative-math \
    -fno-reciprocal-math

ROOT_DIR := $(shell dirname $(PWD))

# Include folders
INCLUDE = -I$(ROOT_DIR)/include -I$(ROOT_DIR)/test

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -L$(ROOT_DIR) -L$(ROOT_DIR)/test -lxatu -larmadillo -lopenblas -llapack -larpack -fopenmp -lgfortran

dummy_bin_folder := $(shell mkdir -p bin)

tests: tests.cpp
	$(CC) -o bin/$@.x $< $(CFLAGS) $(TEST_FP_FLAGS) $(INCLUDE) $(LIBS) 

clean:
	rm -f bin/*.x
