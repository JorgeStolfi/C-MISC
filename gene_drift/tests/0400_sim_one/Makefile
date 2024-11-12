# Last edited on 2023-06-01 07:30:09 by stolfi

PROGDIR := ../..
PROG := gdr_test_sim

OUT_DIR := out

N_TESTS := 1
N_YEARS := 500
INI_SIZE := 50
FIN_SIZE := 50
PARMS := 25 0.920 3 15 60

all: ${PROGDIR}/${PROG} Makefile 
	mkdir -p ${OUT_DIR}
	${PROGDIR}/${PROG} --help
	${PROGDIR}/${PROG} --info
	${PROGDIR}/${PROG} \
          -nTests ${N_TESTS} \
          -nYears ${N_YEARS} \
          -iniSize ${INI_SIZE} \
          -finSize ${FIN_SIZE} \
          -demoParms ${PARMS}
        

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
