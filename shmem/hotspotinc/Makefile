
OSHMEM_CC=cc
OSHMEM_C_FLAGS=-O3 -g
LIBS=

all: hotspotinc

hotspotinc: hotspotinc.c
	$(OSHMEM_CC) $(OSHMEM_C_FLAGS) -o $@ $< $(LIBS)

clean:
	rm hotspotinc
