# Last edited on 2012-12-20 19:00:12 by stolfilocal

PROG := logsolve
PROGDIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

run: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -multiplier 16 16 -factors 21023 51019 \
          -heuristic 2
        
