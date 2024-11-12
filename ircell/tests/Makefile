# Last edited on 2017-03-04 20:20:20 by jstolfi

PROGDIR := ..

PROG := ircell

all: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} > out/cell.ps
	ps2pdf out/cell.ps out/cell.pdf
	evince out/cell.pdf

include ../../../GENERIC-PROGS-TEST.make
