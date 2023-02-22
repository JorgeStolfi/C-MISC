# Last edited on 2023-02-21 20:00:53 by stolfi

PROG := shockwave
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: clean test-run

clean::
	-/bin/rm -f *.eps *.ps

PSVIEW := evince
SHOWPGM := display

OUTNAME := out/sproing

EPSFILES := \
  ${OUTNAME}.eps
  
test-run: ${EPSFILES}

${EPSFILES}: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG}
	for psf in ${EPSFILES}; do \
          ${PSVIEW} $$psf; \
        done
