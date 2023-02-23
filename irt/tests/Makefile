# Last edited on 2023-02-22 20:46:41 by stolfi

PROG := irt_main
PROGDIR := ..

TESTS := \
  torus/400/ia:id:mx:aa
  
NONTESTS := \
  ball/128/ia:id:mx:aa \
  bitorus/200/ia:id:mx:aa \
  cube/200/ia:id:mx:aa \
  equipot/200/ia:id:mx:aa \
  equipot2/200/ia:id:mx:aa \
  equipot3/200/ia:id:mx:aa \
  humps/200/ia:id:mx:aa \
  linktori/200/ia:id:mx:aa \
  linktori2/200/ia:id:mx:aa \
  max/200/ia:id:mx:aa \
  mul/200/ia:id:mx:aa \
  ncubes/200/ia:id:mx:aa \
  plane/200/ia:id:mx:aa \
  plane1/200/ia:id:mx:aa \
  plane2/200/ia:id:mx:aa \
  sqr/200/ia:id:mx:aa \
  torus/200/ia:id:mx:aa \
  bitorus/200/ia:id:mx:aa
  
include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

.PHONY:: \
  run-all clean-all \
  run-shape clean-shape \
  run-arith clean-arith \

SHOWPPM := display
  
all: clean-all run-all

clean-all:
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS=$$t clean-shape ; \
        done

run-all: ${PROGDIR}/${PROG}
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS=$$t run-shape ; \
        done

######################################################################
# Recursive make - caller must define ${TESTARGS}

TESTARGS := TESTARGS.IS.UNDEFINED
ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"

ARGS := ${subst /, ,${TESTARGS}}
TNAME   := ${word 1,${ARGS}}
IMGSIZE := ${word 2,${ARGS}}
ARITHS  := ${subst :, ,${word 3,${ARGS}}}

SUBDIR := ${TNAME}

EXPR_FILE := ${SUBDIR}/${TNAME}.expr
PCODE_FILE := ${SUBDIR}/${TNAME}.pcode
PARMS_FILE := ${SUBDIR}/${TNAME}.parms

clean-shape:
	for m in ${ARITHS} ; do \
          ${MAKE} TESTARGS=${TESTARGS} ARITH=$$m clean-arith ; \
        done

run-shape:
	for m in ${ARITHS} ; do \
          ${MAKE} TESTARGS=${TESTARGS} ARITH=$$m run-arith ; \
        done

OUT_DIR := ${SUBDIR}/out

${OUT_DIR}:
	mkdir -pv ${OUT_DIR}

EXCOMPILER := ${STOLFIHOME}/bin/expr_to_pcode.sh

${PCODE_FILE}: ${EXPR_FILE} ${EXCOMPILER} Makefile
	${EXCOMPILER} < ${EXPR_FILE} > ${PCODE_FILE}

######################################################################
# Recursive make - caller must define ${ARITH}

ifneq "/${ARITH}" "/"

DATE := ${shell date '+%Y-%m-%d-%H%M%S'}

RUN :=  ${TNAME}-${IMGSIZE}-${ARITH}

OUT_PREFIX := ${OUT_DIR}/${RUN}

PPMFILE := ${OUT_PREFIX}.ppm
EVPGMFILE := ${OUT_PREFIX}-evals.pgm
LOGFILE := ${OUT_PREFIX}-${DATE}.log

PRINT_RAY := -1 -1
PLOT_RAY := -1 -1

clean-arith: 
	rm -f ${OUT_DIR}/*.{pgm,ppm,log,eps}

run-arith: ${PPMFILE}
	@echo '= = = did ${PPMFILE} = = = = = = = = = = = = = = = = = = = = ='

${PPMFILE}: ${PCODE_FILE} ${PARMS_FILE} ${PROGDIR}/${PROG} ${OUT_DIR} Makefile
	${PROGDIR}/${PROG} \
            -imageSize ${IMGSIZE} ${IMGSIZE} \
            -pcodeFile ${PCODE_FILE} \
            -arith ${ARITH} \
            -parmsFile ${PARMS_FILE} \
            -printRay ${PRINT_RAY} \
            -plotRay ${PLOT_RAY} \
            -outPrefix ${OUT_PREFIX} \
          2>&1 \
          | tee ${LOGFILE}
	display ${PPMFILE} ${EVPGMFILE}

endif
# end make ${ARITH}
######################################################################

endif
# end make ${TESTARGS}
######################################################################
