
OSHMEM_CC=cc
OSHMEM_C_FLAGS=-O3 -g
LIBS=

all: randominc

randominc: randominc.c
	$(OSHMEM_CC) $(OSHMEM_C_FLAGS) -o $@ $< $(LIBS)

clean:
	rm randominc
