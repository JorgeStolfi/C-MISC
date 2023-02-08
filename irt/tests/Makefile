# Last edited on 2022-10-23 17:50:22 by stolfi

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

SHOWPPM := display
  
all: ${PROGDIR}/${PROG}
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS=$$t single ; \
        done

######################################################################
# Recursive make - caller must define ${TESTARGS}

TESTARGS := TESTARGS.IS.UNDEFINED
ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"

ARGS := ${subst /, ,${TESTARGS}}
TNAME   := ${word 1,${ARGS}}
IMGSIZE := ${word 2,${ARGS}}
ARITHS  := ${word 3,${ARGS}}

SUBDIR := ${TNAME}

EXPR_FILE := ${SUBDIR}/${TNAME}.expr
PCODE_FILE := ${SUBDIR}/${TNAME}.pcode

single: ${PCODE_FILE}
	for m in ${subst :, ,${ARITHS}} ; do \
          ${MAKE} TESTARGS=${TESTARGS} ARITH=$$m single-arith ; \
        done

# ${PCODE_FILE}: ${EXPR_FILE} ${EXCOMPILER}
# 	${EXCOMPILER} < ${EXPR_FILE} > ${PCODE_FILE}

######################################################################
# Recursive make - caller must define ${ARITH}

ARITH := ARITH.IS.UNDEFINED
ifneq "/${ARITH}" "/ARITH.IS.UNDEFINED"

DATE := ${shell date '+%Y-%m-%d-%H%M%S'}

RUN :=  ${TNAME}-${IMGSIZE}-${ARITH}

PARMFILE := ${SUBDIR}/${RUN}.parms
PPMFILE := ${SUBDIR}/out/${RUN}.ppm
EVPGMFILE := ${SUBDIR}/out/${RUN}-evals.pgm
LOGFILE := ${SUBDIR}/out/${RUN}-${DATE}.log

single-arith: ${PPMFILE}
	@echo '= = = did ${PPMFILE} = = = = = = = = = = = = = = = = = = = = ='

${PPMFILE}: ${PCODE_FILE} ${PARMFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} ${SUBDIR} ${RUN} 2>&1 | tee ${LOGFILE}
	display ${PPMFILE} ${EVPGMFILE}

endif
# end make ${ARITH}
######################################################################

endif
# end make ${TESTARGS}
######################################################################
