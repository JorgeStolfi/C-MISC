# Last edited on 2023-06-01 07:30:15 by stolfi

PROGDIR := ../..
PROG := gdr_test_gen_children

OUT_DIR := out

N_TESTS := 20000
PARMS := 25 0.920 3 15 60

all: ${PROGDIR}/${PROG} Makefile 
	mkdir -p ${OUT_DIR}
	${PROGDIR}/${PROG} --help
	${PROGDIR}/${PROG} --info
	${PROGDIR}/${PROG} \
          -nTests ${N_TESTS} \
          -demoParms ${PARMS}
        

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
