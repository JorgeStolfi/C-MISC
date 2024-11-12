# Last edited on 2011-06-09 13:05:47 by stolfi

PROG := tcoptim
PROGDIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

run: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} -size 3 3 -outName out/foo
	${PROGDIR}/${PROG} -size 2 2 -outName out/foo
	${PROGDIR}/${PROG} -size 2 3 -outName out/foo
	${PROGDIR}/${PROG} -size 2 4 -outName out/foo

not-run:
	${PROGDIR}/${PROG} -size 3 4 -maxBifuns 35000000 -maxGates 5 -outName out/foo
