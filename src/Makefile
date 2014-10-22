# Profiling options
#PROFILING=-g -pg

# Compiler warnings
WARN=-Wall

CC = gcc
CFLAGS=-g -std=c99 -O3 -mtune=core2 $(WARN) $(PROFILING)
LINKFLAGS=-g
#LIBS=-lpthread

AR = ar
ARFLAGS = rvs

PROG=salt
VERSION=0.0.1
SLIB = $(PROG).a.$(VERSION)

DEPS=salt.h Makefile

OBJS=query.o util.o maps.o popcount.o overlap_nuc.o \
overlap_nuc4_sse_8.o overlap_nuc4_sse_16.o overlap_nuc4_avx2_8.o \
overlap_nuc4_avx2_16.o

.SUFFIXES:.o .c

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -mavx2 -c -o $@ $<

all: $(SLIB)

$(SLIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *~ $(PROG) gmon.out
