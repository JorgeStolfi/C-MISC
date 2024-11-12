# Last edited on 2023-02-04 07:02:45 by stolfi

SOLVEPROG := SPTimeOnlySolveG3C1

PROG := ${SOLVEPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 32

# mass:friction:spring:A0:A1:tStep:RHS:solution

RUNS := \
  +1.00:00.00:00.00:+1.000:-1.000:1.000:TFself:exnt

NONRUNS := \
  +1.00:00.00:00.00:+1.000:-1.000:1.000:TFexnt:exnt \
  +0.50:+1.00:00.00:+1.000:-1.000:1.000:TFnhff:exnt \
  00.00:+1.00:+1.00:+1.000:-1.000:1.000:TFzero:exnt \
  00.00:+1.00:00.00:+1.000:-1.000:1.000:TFnegf:exnt \
  +1.00:00.00:00.00:+1.000:-1.000:1.000:TFself:exnt \
  +1.00:00.00:+1.00:+1.000:00.000:1.000:TFzero:cost \
  +1.00:+2.00:+2.00:+1.000:-1.000:1.000:TFzero:amor \
  +1.00:00.00:00.00:00.000:00.000:1.000:TFstep:htsq

RHS_RUNS := \
  +1.00:00.00:-1.00:+1.000:-1.000:0.100:TFzero:exnt

TSTEP_RUNS := \
  +1.00:00.00:+1.00:+1.000:00.000:1.000:TFzero:cost \
  +1.00:00.00:+1.00:+1.000:00.000:0.500:TFzero:cost \
   \
  +1.00:00.00:00.00:00.000:00.000:1.000:TFstep:htsq \
  +1.00:00.00:00.00:00.000:00.000:0.100:TFstep:htsq

.PHONY:: run clean-single run-single \
  make-approx show-ps

all: clean run
# all: run

clean::
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" clean-single ; \
        done

run: 
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" run-single ; \
        done
  
# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS  := ${subst :, ,${RUNARGS}}
MASS     := ${word 1,${FIELDS}}
FRICTION := ${word 2,${FIELDS}}
SPRING   := ${word 3,${FIELDS}}
A0       := ${word 4,${FIELDS}}
A1       := ${word 5,${FIELDS}}
TSTEP    := ${word 6,${FIELDS}}
RHSNAME  := ${word 7,${FIELDS}}
SOLNAME  := ${word 8,${FIELDS}}

NSTEPS := 20
TMIN := 0.00000
TMAX := ${shell echo '${TMIN} + ${NSTEPS}*${TSTEP}' | bc -lq}

OUTDIR  := ${DATADIR}/TOSOL/${RHSNAME}
OUTNAME := ${OUTDIR}/s

clean-single: 
	-/bin/rm -f ${OUTNAME}.erp \
          ${OUTNAME}.ps 

ININAME := /tmp/test
INIFRAME := ${ININAME}-000000.tst
  
${INIFRAME}: Makefile
	echo 'begin SPTimeOnlyFunctionFrame (format of 2005-08-18)' > $@
	echo 'time = 0.0000' >> $@
	echo 'coeffs = 2' >> $@
	echo "${A0}" >> $@
	echo "${A1}" >> $@
	echo 'end SPTimeOnlyFunctionFrame' >> $@

run-single: ${INIFRAME} ${OUTDIR} ${PROGDIR}/${SOLVEPROG} Makefile
	${PROGDIR}/${SOLVEPROG} \
            -rhsName ${RHSNAME} \
            -outName ${OUTNAME} \
            -cMass ${MASS} -cFriction ${FRICTION} -cSpring ${SPRING} \
            -iniName ${ININAME} -solName ${SOLNAME} \
            -timeIni ${TMIN} -timeStep ${TSTEP} -nSteps ${NSTEPS} \
            -plotSteps 4 \
            -maxIter 20 -absTol 0.0000000001 -relTol 0.0000000001 \
            -smpOrder ${SMPORDER} -verbose
	show-plot -show ${OUTNAME} ${TMIN} ${TMAX}

${OUTDIR}:
	mkdir -p ${OUTDIR}

# PSVIEW := kghostview
PSVIEW := evince

show-ps: ${OUTNAME}.ps
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-{app,err}-{f,v}.eps` ; do \
          ${PSVIEW} $$el ; \
        done

endif
