
CC=mpicc
CFLAGS=-O3 -std=c99
EXTRA_CFLAGS=-D_POSIX_C_SOURCE=200809L
LIBS=

sweep3d: sweep3d.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -o $@ $< $(LIBS)

clean:
	rm sweep3d

