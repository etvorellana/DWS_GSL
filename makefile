
CC = gcc
CFLAGS = -O2 -Wall

all: dws

libdws: libdws.c dws.h
	$(CC) $(CFLAGS) -c libdws.c

dws: DWS_C.c dws.h libdws.o
	$(CC) $(CFLAGS) -o dws -lm -lgsl -lgslcblas DWS_C.c libdws.o

clean:
	rm libdws.o dws
