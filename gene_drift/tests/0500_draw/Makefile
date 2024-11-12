# Last edited on 2023-06-02 09:01:22 by stolfi

PROGDIR := ../..
PROG := gdr_test_draw

DEBUG := BIG

Y_START := +1750 # First year.
Y_STOP :=  +1950 # Last year.
TESTS :=   2.A # Cohort sizes and file tags.

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/gdr

all: 
	mkdir -p ${OUT_DIR}
	rm -fv ${OUT_PREFIX}-*.eps
	for tst in ${TESTS} ; do \
	  ${MAKE} TEST=$${tst}  single ; \
	done

ifneq "/${TEST}" "/"
# ######################################################################
# Recursive section ${TEST}

ARGS := ${subst ., , ${TEST}}
COHORT_SIZE := ${word 1, ${ARGS}}
TEX_TAG := ${word 2, ${ARGS}}

# Cohort sizes:
INI_SIZE := ${COHORT_SIZE}
FIN_SIZE := ${COHORT_SIZE} 

OUT_FULL_PREFIX := ${OUT_PREFIX}-${TEX_TAG}

OUT_PLOT_FILE := ${OUT_FULL_PREFIX}-0-01-evol.eps

single: ${OUT_PLOT_FILE} 
	evince ${OUT_FULL_PREFIX}-?-??-evol.eps

CH_PARMS_M := 15 0.810 3 10 30
CH_PARMS_F :=  6 0.180 3 10 30

${OUT_PLOT_FILE}: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          -iniSize ${INI_SIZE} \
          -finSize ${FIN_SIZE} \
          -yStart ${Y_START} \
          -yStop ${Y_STOP} \
          -demoParms 0 ${CH_PARMS_M} \
          -demoParms 1 ${CH_PARMS_F} \
          -outPrefix ${OUT_PREFIX} \
          -tag ${TEX_TAG}

# End section ${COHORT_SIZE}
# ######################################################################
endif


include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
