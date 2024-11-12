# Last edited on 2023-06-03 15:28:34 by stolfi

PROGDIR := ../..
PROG := gdr_main

DEBUG := BIG

ifeq "/${DEBUG}" "/TINY"
  Y_START := +1750 # First year.
  Y_STOP :=  +1950 # Last year.
  N_RUNS :=      8 # Number of runs per sex.
  TESTS :=   10.A 40.B # Cohort sizes and file tags.
else ifeq "/${DEBUG}" "/SMALL"
  Y_START := +1550 # First year.
  Y_STOP :=  +1950 # Last year.
  N_RUNS :=     50 # Number of runs per sex.
  TESTS :=   20.A 80.B # Cohort sizes and file tags.
else ifeq "/${DEBUG}" "/BIG"
  Y_START := -7950 # First year.
  Y_STOP :=  +1950 # Last year.
  N_RUNS :=   1000 # Number of runs per sex.
  TESTS := 400.B 100.A # Cohort sizes and file tags.
endif

COHORT_SIZES := ${shell echo ${TESTS} | sed -e 's:[.].*$$::g'}

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/gdr

all: 
	for tst in ${TESTS} ; do \
	  ${MAKE} TEST=$${tst}  single ; \
	done
	./plot_ratios.sh ${OUT_PREFIX} ${TESTS}
	./plot_ratios_spreads.sh ${OUT_PREFIX} ${TESTS}

ifneq "/${TEST}" "/"
# ######################################################################
# Recursive section ${TEST}

ARGS := ${subst ., , ${TEST}}
COHORT_SIZE := ${word 1, ${ARGS}}
TEX_TAG := ${word 2, ${ARGS}}

# Cohort sizes:
NI_INI := ${COHORT_SIZE}
NI_FIN := ${COHORT_SIZE} 

OUT_FULL_PREFIX := ${OUT_PREFIX}-${TEX_TAG}

OUT_RATIO_FILE := ${OUT_FULL_PREFIX}-ratio.txt
OUT_PARMS_FILE := ${OUT_FULL_PREFIX}-parms.tex
OUT_PROBS_FILE := ${OUT_FULL_PREFIX}-probs.tex
  
single: ${OUT_RATIO_FILE} plot_nlins.sh
	./plot_nlins.sh       ${OUT_FULL_PREFIX} ${COHORT_SIZE}
	./plot_populations.sh ${OUT_FULL_PREFIX} ${COHORT_SIZE}
	./plot_rsamp.sh       ${OUT_FULL_PREFIX} ${COHORT_SIZE}

CH_PARMS_M := 25 0.810 3 15 65
CH_PARMS_F :=  8 0.180 3 15 45

${OUT_RATIO_FILE}: ${PROGDIR}/${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROGDIR}/${PROG} \
          -nRuns ${N_RUNS} \
          -iniSize ${NI_INI} \
          -finSize ${NI_FIN} \
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
