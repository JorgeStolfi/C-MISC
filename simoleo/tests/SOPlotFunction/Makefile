# Last edited on 2023-02-04 07:04:12 by stolfi

PROG := SOPlotFunction

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

RUNS := \
  teste0:01.001:00.100:10.0:linx:ps: \
  teste1:01.001:00.100:01.0:wave-p004-m003:ps:-treeName%tree/circle-min%-noBands% \
  teste2:02.001:00.100:01.0:tent1:ps:-noIsolines% \
  teste3:00.251:00.025:01.0:mhat:eps:% \
  teste5:01.001:00.100:01.0:gaus:eps:-noIsolines% \
  teste6:02.001:00.100:01.0:dgdx:eps:%

RUNS_LATER := \
  unit \
  linx \
  liny \
  lind \
  sqrx \
  para \
  hump \
  sinx \
  sxsy \
  sxcy \
  cbsm

NON_RUNS := \
  linz

.PHONY:: run \
  clean-single run-single

all: run

clean:: 
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" clean-single ; \
        done
  
run: 
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" run-single ; \
        done

######################################################################
# For recursive "make" only -- caller must define ${RUNARGS}
# 

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS    := ${subst :, ,${RUNARGS}}
OUTNAME   := ${word 1,${FIELDS}}
FRANGE    := ${word 2,${FIELDS}}
FSTEP     := ${word 3,${FIELDS}}
MESHSIZE  := ${word 4,${FIELDS}}
FUNCNAME  := ${word 5,${FIELDS}}
OUTEXT    := ${word 6,${FIELDS}}
OPTSX     := ${word 7,${FIELDS}}

OUTPREFIX :=  out/${OUTEXT}/${OUTNAME}

OPTS := ${subst %, ,${OPTSX}}

FUNDIR := fun
PSNAME := ${OUTPREFIX}.${OUTEXT}

PSVIEW := evince

run-single: ${PSNAME}

${PSNAME}: ${FUNDIR}/${FUNCNAME}.fun ${PROGDIR}/SOPlotFunction
	${PROGDIR}/SOPlotFunction \
          -funcName ${FUNDIR}/${FUNCNAME} \
          -fRange ${FRANGE} -fStep ${FSTEP} \
          -${OUTEXT} ${OPTS} \
          -meshSize ${MESHSIZE} \
          -outName ${OUTPREFIX}
	${PSVIEW} ${PSNAME}

endif
# end ${RUNARGS}
######################################################################
