# Last edited on 2008-06-10 20:44:08 by stolfi

MAKEPROG := SPComputeErrorMap
PLOTPROG := SPPlotFunction

PROGS := \
  ${MAKEPROG} \
  ${PLOTPROG}
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 40
SO := ${SMPORDER}

SHOW := 1

# ABASIS:HARMDEG:HARMEXC:HOMOFLAG

RUNS := \
  pw/icosa/raw/nh-d6c1:17:06:@ \
  pw/icosa/raw/h-d7c1:17:07:-homogeneous

RUNS_TO_DO :=

.PHONY:: run clean-single run-single

# all: clean run
# all: run
all: run

clean:
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" clean-single ; \
        done

run:
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" run-single ; \
        done

# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS    := ${subst :, ,${RUNARGS}}
ABASIS    := ${word 1,${FIELDS}}
HARMDEG   := ${word 2,${FIELDS}}
HARMEXCL  := ${word 3,${FIELDS}}
HOMOFLAG  := ${subst @,,${word 4,${FIELDS}}}

ABASNAME := ${DATADIR}/sp/${ABASIS}/b
AMATDIR := ${DATADIR}/sp/${ABASIS}
AMATNAME := ${AMATDIR}/b

OUTNAME := ${ABASNAME}-emap

clean-single:
	-/bin/rm -f ${OUTNAME}{.sfn,.ps,-??????.eps}
	-/bin/rm -f ${OUTNAME}.log

run-single: make-error-map plot-error-map

make-error-map: ${OUTNAME}.sfn

plot-error-map: ${OUTNAME}.ps

${AMATDIR}:
	mkdir -p ${AMATDIR}


${OUTNAME}.mat: ${PROGDIR}/${MAKEPROG} Makefile \
            ${AMATDIR} ${ABASNAME}.bas
	time ${PROGDIR}/${MAKEPROG} \
            -basisName ${ABASNAME} \
            -matName ${AMATNAME} \
            -outName ${OUTNAME} \
            -smpOrder ${SMPORDER} \
            -degree ${HARMDEG} -exclude ${HARMEXCL} \
            ${HOMOFLAG} \
          2>&1 | tee ${OUTNAME}.log
	${PSVIEW} ${OUTNAME}.ps

PSOPTIONS := \
  -ps -paperSize letter \
  -figSize 80 -meshSize 0.666 \
  -caption "${ABASIS}:${HARMDEG}:${HARMEXCL}" \
  -caption "step = %S"

${OUTNAME}.ps: ${PROGDIR}/${PLOTPROG} Makefile \
            ${OUTDIR} ${OUTNAME}.sfn
	time ${PROGDIR}/${PLOTPROG} \
            -funcName ${OUTNAME} \
            -outName ${OUTNAME} \
            -showTriang -bothSides \
	    -${PSOPTIONS} \
	    -obs 1 14 09 07 -radius 1.2 \
            -autoRange 0.9 \
          2>&1 | tee ${OUTNAME}.log
	if [ ${SHOW} ] ; then \
	  ${PSVIEW} ${OUTNAME}.ps ; \
	fi

endif
