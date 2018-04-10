BINDIR=bin
CC=gcc
CFLAGS=-Wall -Wextra -pedantic -O2

.PHONY: all
all: $(BINDIR)/dice-opt $(BINDIR)/dice-opt-omp

$(BINDIR)/dice-opt: dice-opt.c
	$(CC) $(CFLAGS) -Wno-unknown-pragmas -o $@ $^

$(BINDIR)/dice-opt-omp: dice-opt.c
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

.PHONY: clean
clean:
	rm $(BINDIR)/dice-opt $(BINDIR)/dice-opt-omp
