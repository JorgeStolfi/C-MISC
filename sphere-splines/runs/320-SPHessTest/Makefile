# Last edited on 2023-02-04 07:03:16 by stolfi

MAKEPROG := SPHessTest

PROGS := ${MAKEPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data
ALLFUNDIR := ${DATADIR}/fn

FUNS := \
  pw/h-d7c1-00001
  
FUNS_BUGGY := \
  proc/hr75
  
FUNS_LATER := \
  proc/unit \
  proc/linx \
  proc/liny \
  proc/linz \
  proc/sqrx \
  proc/sqry \
  proc/cubx \
  proc/qrtx \
  proc/qrtz \
  proc/quix \
  proc/sexx \
  proc/sepx \
  proc/octx \
  proc/band \
  proc/hrm2 \
  proc/hr75 \
  proc/loca \
  proc/fxyz \
  proc/expx \
  proc/expz \
  proc/sinx \
  proc/cosx \
  proc/lrry \
  proc/rain \
  proc/spir

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
FUNNAME := ${FUNDIR}/f
FUNFILE := ${FUNNAME}.sfn
OUTNAME := ${FUNNAME}-hess-error

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
  -caption "${FUN} [%S]" \
  -caption "step = %D"

SHOWPSFILES := ${PSVIEW} ${PSFILE}

endif
#
######################################################################

single: ${PSTARGET}

${PSTARGET}: \
            ${FUNFILE} \
            ${PROGDIR}/SPHessTest \
            Makefile
	${PROGDIR}/SPHessTest \
          -sfnName ${FUNNAME} \
          -smpOrder 10 \
          -delta 0.000001 \
          -outName ${OUTNAME} \
          ${PSOPTIONS} \
          -obs 1 14 09 07 -radius 1.2 \
          -autoRange 0.9 -fRange 1.0 -fStep 0.2
	if [ ${SHOW} -gt 0 ]; then ${SHOWPSFILES}; fi


endif
# ${FUN}
######################################################################
