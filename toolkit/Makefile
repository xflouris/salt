# Profiling options
#PROFILING=-g -pg

# Compiler warnings
WARN=-Wall

CC = gcc
INCDIR = ../src
LIBDIR = ../src
CFLAGS=-g -std=c99 -O3 -mtune=core2 -I $(INCDIR) -L $(LIBDIR) $(WARN) $(PROFILING)
LINKFLAGS=-g
LIBS=-lsalt

PROG=salt

DEPS = salt.h Makefile

OBJS = salt.o

.SUFFIXES:.o .c

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -mavx2 -c -o $@ $<

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LIBS)

clean:
	rm -f *.o *~ $(PROG) gmon.out
