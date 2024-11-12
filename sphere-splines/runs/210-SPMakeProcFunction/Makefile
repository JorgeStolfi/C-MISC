# Last edited on 2023-02-04 07:03:45 by stolfi

MAKEPROG := SPMakeProcFunction
PLOTPROG := SPPlotFunction

PROGS := ${MAKEPROG} ${PLOTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data
ALLFUNDIR := ${DATADIR}/fn/proc

FUNS := \
  unit \
  linx \
  liny \
  linz \
  sqrx \
  sqry \
  cubx \
  qrtx \
  qrtz \
  quix \
  sexx \
  sepx \
  octx \
  band \
  hrm2 \
  loca \
  fxyz \
  expx \
  expz \
  sinx \
  cosx \
  lrry \
  rain \
  spir \
  hr75
  
FUNS_LATER :=

SHOW := 1
FORMAT := ps

.PHONY:: single

all: 
	for f in ${FUNS}; do \
          ${MAKE} FUN=$$f single; \
        done

# PSVIEW := kghostview
PSVIEW := evince

######################################################################
# Recursive make - caller must define ${FUN}

FUN := FUN.IS.UNDEFINED
ifneq "/${FUN}" "/FUN.IS.UNDEFINED"

FUNDIR := ${ALLFUNDIR}/${FUN}
OUTNAME := ${FUNDIR}/f
FUNFILE := ${OUTNAME}.sfn

######################################################################
# ${FORMAT} == "eps"

ifeq "/${FORMAT}" "/eps"

EPSFILE_F := ${OUTNAME}-f-000001.eps
EPSFILE_R := ${OUTNAME}-r-000001.eps

PSTARGET := ${EPSFILE_F}

PSOPTIONS := \
  -eps \
  -figSize 150 -meshSize 1.5

SHOWPSFILES := ${PSVIEW} ${EPSFILE_F} ; ${PSVIEW} ${EPSFILE_R}

endif
#
######################################################################

######################################################################
# ${FORMAT} == "ps"

ifeq "/${FORMAT}" "/ps"

PSFILE := ${OUTNAME}.ps
PSFILES := ${PSFILE}

PSTARGET := ${PSFILE}

PSOPTIONS := \
  -ps -paperSize letter \
  -figSize 80 -meshSize 0.666 \
  -caption "${FUN} [%S]"

SHOWPSFILES := ${PSVIEW} ${PSFILE}

endif
#
######################################################################

single: ${FUNFILE} ${PSTARGET}

${FUNFILE}: ${PROGDIR}/${MAKEPROG}
	mkdir -p ${FUNDIR}
	${PROGDIR}/${MAKEPROG} \
          -function ${FUN} \
          -outName ${OUTNAME}

${PSTARGET}: \
            ${FUNFILE} \
            ${PROGDIR}/${PLOTPROG} \
            Makefile
	${PROGDIR}/${PLOTPROG} \
          -funcName ${OUTNAME} \
          -outName ${OUTNAME} \
          ${PSOPTIONS} \
          -obs 1 14 09 07 -radius 1.2 -bothSides \
          -autoRange 0.9 -fRange 1.0 -fStep 0.2
	if [ ${SHOW} -gt 0 ]; then ${SHOWPSFILES}; fi

endif
# ${FUN}
######################################################################
