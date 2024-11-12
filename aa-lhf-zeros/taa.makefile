# Last edited on 2009-02-09 21:09:24 by stolfi

# NONSTANDARD MAKEEFILE

PROGS := taa

N=9
N=20
N=32

CC=/usr/bin/gcc
CFLAGS := -O2 -Wall -DAA_N=$N -g -ffloat-store -frounding-math

GP=gp.o dv.o

LIBS= ${GP} -L/usr/X11R6/lib -lX11 -lm

all:	build install

build: ${PROGS}

install:

uninstall:

clean::
	rm -f *.o aa core c2 a.out

run:	taa
	taa

taa:	${PROG}.o aasolve.h aasolve.o fcircle.c f2circ.c aa.o jr.o ${GP}
	${CC} ${CFLAGS} ${IFLAGS} -o $@ ${PROG}.o aa.o jr.o aasolve.o ${LIBS}

${PROG}.o:	${PROG}.c aasolve.h fcircle.c f2circ.c fcirhex.c jr.h aa.h
	${CC} ${CFLAGS} ${IFLAGS} -o $@ -c ${PROG}.c

aasolve.o: aasolve.c aasolve.c aa.h
	${CC} ${CFLAGS} ${IFLAGS} -o $@ -c aasolve.c

gp:
	${CC} -O2 -c gp.c dv.c
