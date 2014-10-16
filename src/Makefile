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

OBJS=salt.o query.o util.o maps.o score.o popcount.o gen_test.o \
overlap_plain.o overlap_sse_8bit.o overlap_sse_16bit.o overlap_avx2_8bit.o overlap_avx2_16bit.o

.SUFFIXES:.o .c

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -mavx2 -c -o $@ $<

all: $(PROG)
	@echo "Done."

salt: $(OBJS)
	$(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *~ $(PROG) gmon.out
