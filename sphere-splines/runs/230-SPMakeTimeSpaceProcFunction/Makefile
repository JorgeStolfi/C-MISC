# Last edited on 2023-02-04 07:03:30 by stolfi

MAKEPROG := SPMakeTimeSpaceProcFunction
PLOTPROG := SPPlotTimeSpaceFunction

PROGS := ${MAKEPROG} ${PLOTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data
ALLFUNDIR := ${DATADIR}/tsfn/proc

FUNS_VARYING := \
  turn \
  ripl \
  spir \
  rain \
  lint \
  sqtm \
  cost \
  xent \
  xsnt
  
FUNS_CONSTANT := \
  unit \
  linx \
  liny \
  linz \
  sqrx \
  sqry \
  hent

FORMAT := ps

.PHONY:: all-var all-const single

all: all-var all-const

all-var:
	for f in ${FUNS_VARYING}; do \
          ${MAKE} FUN=$$f PLOT=1 SHOW=1 single; \
        done

all-const:
	for f in ${FUNS_CONSTANT}; do \
          ${MAKE} FUN=$$f PLOT=0 SHOW=0 single; \
        done

# PSVIEW := kghostview
PSVIEW := evince

######################################################################
# Recursive make - caller must define ${FUN}

FUN := FUN.IS.UNDEFINED
ifneq "/${FUN}" "/FUN.IS.UNDEFINED"

FUNDIR := ${ALLFUNDIR}/${FUN}
OUTNAME := ${FUNDIR}/f
FUNFILE := ${OUTNAME}.tsf

${FUNFILE}: \
            ${PROGDIR}/${MAKEPROG}
	mkdir -p ${FUNDIR}
	${PROGDIR}/${MAKEPROG} \
          -function ${FUN} \
          -outName ${OUTNAME}

######################################################################
# ${PLOT} == 1

ifeq "/${PLOT}" "/1"

######################################################################
# ${FORMAT} == "eps"

ifeq "/${FORMAT}" "/eps"

EPSFILE_F := ${OUTNAME}-f-000001.eps
EPSFILE_R := ${OUTNAME}-r-000001.eps

PSTARGET := ${EPSFILE_F}

TIMEOPTIONS := \
  -timeStart 0.0000 -timeStep 0.1000 -nTimes 11

PSOPTIONS := \
  -eps \
  -figSize 150 -meshSize 1.5 \
  -bothSides

SHOWPSFILES := ${PSVIEW} ${EPSFILE_F} ; ${PSVIEW} ${EPSFILE_R}

endif
# ${FORMAT}
######################################################################

######################################################################
# ${FORMAT} == "ps"

ifeq "/${FORMAT}" "/ps"

PSFILE := ${OUTNAME}.ps

PSTARGET := ${PSFILE}

TIMEOPTIONS := \
  -timeStart 0.0000 -timeStep 0.1000 -nTimes 35

PSOPTIONS := \
  -ps -paperSize a3 \
  -figSize 43 -meshSize 0.666 \
  -caption "t = %T"

SHOWPSFILES := ${PSVIEW} ${PSFILE}

endif
# ${FORMAT}
######################################################################

${PSTARGET}: \
            ${FUNFILE} \
            ${PROGDIR}/${PLOTPROG} \
            Makefile
	${PROGDIR}/${PLOTPROG} \
          -funcName ${OUTNAME} \
          -outName ${OUTNAME} \
          ${PSOPTIONS} \
          ${TIMEOPTIONS} \
          -obs 1 14 09 07 -radius 1.2 \
          -autoRange 0.9 -fRange 1.0 -fStep 0.2
	if [ ${SHOW} -gt 0 ]; then ${SHOWPSFILES}; fi
 
single: ${FUNFILE} ${PSTARGET}

endif
# ${PLOT}
######################################################################

######################################################################
# ${PLOT} == 0

ifeq "/${PLOT}" "/0"

single: ${FUNFILE}

endif
# ${PLOT}
######################################################################

endif
# ${FUN}
######################################################################
