# Last edited on 2016-09-27 12:13:27 by stolfilocal

PROG := dynstolfi
PROG_DIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make


run: ${PROG_DIR}/${PROG}
	${PROG_DIR}/${PROG} GL_Phi_u 50   0.0 0.5 0.2  1.2 -1.2 2.0 -6.0
        
clean::

actions:
	${MAKE} ${ACTIONS}
        
