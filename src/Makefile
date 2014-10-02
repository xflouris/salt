# Profiling options
#PROFILING=-g -pg

# Compiler warnings
WARN=-Wall

CC = gcc
CFLAGS=-std=c99 -O3 -msse4.1 -mtune=core2 $(WARN) $(PROFILING)
LINKFLAGS=-g
#LIBS=-lpthread

PROG=salt

DEPS=salt.h Makefile

OBJS=salt.o query.o util.o maps.o score.o popcount.o

.SUFFIXES:.o .cc

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROG)

salt: $(OBJS)
	$(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *~ $(PROG) gmon.out
