BINDIR=bin
CC=gcc
CFLAGS=-Wall -Wextra -pedantic -O2
TARGETS= \
	$(BINDIR)/dice-opt \
	$(BINDIR)/dice-opt-omp \
	$(BINDIR)/dice-opt-2 \
	$(BINDIR)/dice-opt-2-omp

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

.PHONY: clean
.PHONY: clean
clean:
	rm $(TARGETS)
