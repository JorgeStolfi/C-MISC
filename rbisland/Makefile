# Last edited on 2011-08-04 00:04:26 by stolfilocal

PROG := rbisland

JS_LIBS := 

all: run
  
include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

POP := 100
GEN := 200
FST := 1
TRY := 300

DAT_FILE := /tmp/${PROG}.dat
PLOT_FILE := out/evol-plot.png

run: ${PROG}
	${PROG} ${POP} ${GEN} ${FST} ${TRY} > ${DAT_FILE}
	./plot-evol.sh ${POP} ${GEN} ${DAT_FILE} > ${PLOT_FILE}
        
        
