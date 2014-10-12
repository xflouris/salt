# Profiling options
#PROFILING=-g -pg

# Compiler warnings
WARN=-Wall

CC = gcc
CFLAGS=-g -std=c99 -O3 -mtune=core2 $(WARN) $(PROFILING)
LINKFLAGS=-g
#LIBS=-lpthread

PROG=salt

DEPS=salt.h Makefile

OBJS=salt.o query.o util.o maps.o score.o popcount.o overlap_plain.o \
overlap_plain_vec.o overlap_plain_sse8.o overlap_plain_sse16.o gen_test.o

.SUFFIXES:.o .c

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -mavx2 -c -o $@ $<


all: $(PROG)

salt: $(OBJS)
	$(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *~ $(PROG) gmon.out
