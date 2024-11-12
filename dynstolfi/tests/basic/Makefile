# Last edited on 2012-07-22 14:30:56 by stolfilocal

PROG := dynstolfi
PROG_DIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make


run: ${PROG_DIR}/${PROG}
	${PROG_DIR}/${PROG} sqout    50
	${PROG_DIR}/${PROG} cirout   50
	${PROG_DIR}/${PROG} cirin    50
	${PROG_DIR}/${PROG} contract 50
	${PROG_DIR}/${PROG} sqwhorl  50
	${PROG_DIR}/${PROG} triwhorl 50
	${PROG_DIR}/${PROG} sqturn   50
	${PROG_DIR}/${PROG} rot144   50
	${PROG_DIR}/${PROG} rotphi   50
	${PROG_DIR}/${PROG} trinop   50
	${PROG_DIR}/${PROG} sqnop    50
	${PROG_DIR}/${PROG} haoui    50
	${PROG_DIR}/${PROG} henon    50
        
clean::

actions:
	${MAKE} ${ACTIONS}
        
