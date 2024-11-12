# Last edited on 2023-02-04 07:03:03 by stolfi

PLOTPROG := SPPlotFunction

PROGS := ${PLOTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

PROGDIR := ../../progs
DATADIR := ../../data

# FUNCTION:TRIANG:PROJLONLAT:SHOWGRID

RUNS := \
  pw/icosa/raw/h-d7c1-OLD/i0099:YES:NO:NO \

RUNS_FOR_SLIDES := \
  harm/d006p003:NO:NO:NO \
  harm/d006p003:NO:YES:NO \
  harm/d006p003:NO:YES:YES \
  proc/lrry:NO:NO:NO \
  proc/lrry:NO:NO:YES \
  proc/lrry:NO:YES:YES \
  proc/zero:NO:NO:YES \
  proc/zero:NO:YES:YES \
  pw/icosa/raw/h-d7c1-OLD/i0015:YES:NO:NO \
  pw/icosa/raw/h-d7c1-OLD/i0015:YES:YES:NO \
  pw/icosa/raw/nh-d6c1/lrry:YES:NO:NO \
  pw/icosa/raw/nh-d6c1/lrry:YES:YES:NO \
  pw/icosa/raw/zero:YES:NO:NO \
  pw/icosa/raw/zero:YES:YES:NO

RUNS_PROC := \
  proc/unit:NO:NO:NO \
  proc/linx:NO:NO:NO \
  proc/liny:NO:NO:NO \
  proc/linz:NO:NO:NO \
  proc/sqrx:NO:NO:NO \
  proc/sqry:NO:NO:NO \
  proc/cubx:NO:NO:NO \
  proc/qrtx:NO:NO:NO \
  proc/qrtz:NO:NO:NO \
  proc/quix:NO:NO:NO \
  proc/sexx:NO:NO:NO \
  proc/sepx:NO:NO:NO \
  proc/octx:NO:NO:NO \
  proc/band:NO:NO:NO \
  proc/hrm2:NO:NO:NO \
  proc/loca:NO:NO:NO \
  proc/fxyz:NO:NO:NO \
  proc/expx:NO:NO:NO \
  proc/expz:NO:NO:NO \
  proc/sinx:NO:NO:NO \
  proc/cosx:NO:NO:NO \
  proc/lrry:NO:NO:NO \
  proc/rain:NO:NO:NO \
  proc/spir:NO:NO:NO \
  proc/hr75:NO:NO:NO

RUNS_HARMONIC := \
  harm/d001n001:NO:NO:NO

RUNS_SPLINES := \
  pw/icosa/raw/h-d7c1-OLD/i0000:YES:YES

.PHONY:: single

all: 
	for r in ${RUNS}; do \
          ${MAKE} RUNARGS=$$r single; \
        done

FORMAT := eps
SHOW := 1
# PSVIEW := kghostview
PSVIEW := evince

######################################################################
# Recursive make - caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS := ${subst :, ,${RUNARGS}}
FUNC         := ${word 1,${FIELDS}}
SHOWTRI      := ${word 2,${FIELDS}}
PROJLONLAT   := ${word 3,${FIELDS}}
SHOWGRID     := ${word 4,${FIELDS}}

ALLFUNCDIR := ${DATADIR}/fn

FUNCDIR := ${ALLFUNCDIR}/${FUNC}
FUNCNAME := ${FUNCDIR}/f
FUNCFILE := ${FUNCNAME}.sfn

PROJTAG := ${subst YES,-lon-lat,${subst NO,,${PROJLONLAT}}}
GRIDTAG := ${subst YES,-grd,${subst NO,,${SHOWGRID}}}

OUTNAME := ${FUNCDIR}/f${PROJTAG}${GRIDTAG}

ifeq "/${PROJLONLAT}" "/NO"
  PROJOPTIONS := \
    -obs 1 14 09 07 -radius 1.2 -bothSides \
    -light 3 4 5 -shadow 0.1
else
  PROJOPTIONS := \
    -projLonLat -radius 1.2 -shadow 0
endif

######################################################################
# ${FORMAT} == "eps"

ifeq "/${FORMAT}" "/eps"

PSOPTIONS := \
  -eps \
  -figSize 150 -meshSize 1.5

ifeq "/${PROJLONLAT}" "/NO"
  EPSFILE_F := ${OUTNAME}-fn-f.eps
  EPSFILE_R := ${OUTNAME}-fn-r.eps
  PSTARGET := ${EPSFILE_F}
  SHOWPSFILES := ${PSVIEW} ${EPSFILE_F} ; ${PSVIEW} ${EPSFILE_R}
else
  EPSFILE_F := ${OUTNAME}-fn-f.eps
  PSTARGET := ${EPSFILE_F}
  SHOWPSFILES := ${PSVIEW} ${EPSFILE_F}
endif

endif
#
######################################################################

######################################################################
# ${FORMAT} == "ps"

ifeq "/${FORMAT}" "/ps"

ifeq "/${PROJLONLAT}" "/NO"
  PSOPTIONS := \
    -eps -paperSize letter \
    -figSize 80 -meshSize 0.666 \
    -caption "${FUNC} [%S]" \
    -caption "step = %D"
else
  PSOPTIONS := \
    -eps -paperSize letter \
    -figSize 160 -meshSize 1.000 \
    -caption "${FUNC}" \
    -caption "step = %D"
endif

PSFILE := ${OUTNAME}.ps
PSFILES := ${PSFILE}

PSTARGET := ${PSFILE}
  
SHOWPSFILES := ${PSVIEW} ${PSFILE}

endif
#
######################################################################

single: ${PSTARGET}

${PSTARGET}: \
            ${FUNCNAME}.sfn \
            ${PROGDIR}/${PLOTPROG} \
            Makefile
	${PROGDIR}/${PLOTPROG} \
          -funcName ${FUNCNAME} \
          -outName ${OUTNAME} \
          ${subst YES,-showTriang,${subst NO,,${SHOWTRI}}} \
          ${subst YES,-gridDots 40 20,${subst NO,,${SHOWGRID}}} \
          ${PSOPTIONS} \
          ${PROJOPTIONS} \
          -autoRange 0.9 -fRange 1.0 -fStep 0.2
	if [ ${SHOW} -gt 0 ]; then ${SHOWPSFILES}; fi


endif
# ${RUNARGS}
######################################################################
