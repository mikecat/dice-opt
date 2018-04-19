BINDIR=bin
CC=gcc
CFLAGS=-Wall -Wextra -pedantic -O2
GPUC=nvcc
GPUCFLAGS=
TARGETS= \
	$(BINDIR)/dice-opt \
	$(BINDIR)/dice-opt-omp \
	$(BINDIR)/dice-opt-2 \
	$(BINDIR)/dice-opt-2-omp \
	$(BINDIR)/dice-opt-2-prank \
	$(BINDIR)/dice-opt-2-prank-omp \
	$(BINDIR)/dice-opt-gpu \
	$(BINDIR)/dice-opt-gpu-shared

.PHONY: all
all: $(TARGETS)

$(BINDIR)/dice-opt: dice-opt.c
	$(CC) $(CFLAGS) -Wno-unknown-pragmas -o $@ $^

$(BINDIR)/dice-opt-omp: dice-opt.c
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

$(BINDIR)/dice-opt-2: dice-opt-2.c
	$(CC) $(CFLAGS) -Wno-unknown-pragmas -o $@ $^

$(BINDIR)/dice-opt-2-omp: dice-opt-2.c
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

$(BINDIR)/dice-opt-2-prank: dice-opt-2-prank.c
	$(CC) $(CFLAGS) -Wno-unknown-pragmas -o $@ $^

$(BINDIR)/dice-opt-2-prank-omp: dice-opt-2-prank.c
	$(CC) $(CFLAGS) -fopenmp -o $@ $^

$(BINDIR)/dice-opt-gpu: dice-opt-gpu.cu
	$(GPUC) $(GPUCFLAGS) -o $@ $^

$(BINDIR)/dice-opt-gpu-shared: dice-opt-gpu-shared.cu
	$(GPUC) $(GPUCFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm $(TARGETS)
