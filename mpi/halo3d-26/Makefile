
CC=mpicc
CFLAGS=-O3 -std=c99
EXTRA_CFLAGS=-D_POSIX_C_SOURCE=200809L
LIBS=

halo3d-26: halo3d-26.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -o $@ $< $(LIBS)

clean:
	rm halo3d-26

