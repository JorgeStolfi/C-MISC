# Last edited on 2012-07-22 14:35:24 by stolfilocal

PROG  := unitexture

PROGDIR  := ../..
  
include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

.PHONY:: run

all: clean run show

OUTNAME := out

clean:: 
	-/bin/rm -f ${OUTNAME}-*.*
  
run: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
          -size 8 \
          -iterations 200
