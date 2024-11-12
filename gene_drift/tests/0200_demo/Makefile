# Last edited on 2023-06-01 07:30:27 by stolfi

PROGDIR := ../..
PROG := gdr_test_demo

OUT_DIR := out
OUT_PREFIX := ${OUT_DIR}/gdr-tcp
TEX_TAG := A
OUT_FULL_PREFIX := ${OUT_PREFIX}-${TEX_TAG}

OUT_PROBS_M_FILE := ${OUT_FULL_PREFIX}-0-probs.tex
OUT_PROBS_F_FILE := ${OUT_FULL_PREFIX}-1-probs.tex

OUT_FILES := ${OUT_PROBS_M_FILE} ${OUT_PROBS_F_FILE}

all: ${OUT_PROBS_M_FILE}
  
PARMS_M := 25 0.920 3 15 60
PARMS_F :=  8 0.350 3 15 45

${OUT_PROBS_M_FILE}: ${PROGDIR}/${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROGDIR}/${PROG} --help
	${PROGDIR}/${PROG} --info
	${PROGDIR}/${PROG} \
          -demoParms ${PARMS_M} \
          -demoParms ${PARMS_F}
	for f in ${OUT_FILES} ; do \
          echo "=== $$f ===" 1>&2 ; \
          cat $$f 1>&2 ; \
        done
        

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
