# Last edited on 2023-02-04 07:03:25 by stolfi

MAKEPROG := SPMakeHarmonicFunction
PLOTPROG := SPPlotFunction

PROGS := ${MAKEPROG} ${PLOTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data
ALLFUNDIR := ${DATADIR}/fn/harm

# FUNCTAG:DEGREE:ORDER

RUNS := \
  d000o000:000:0000 \
  d001n001:001:-001 \
  d001o000:001:0000 \
  d001p001:001:+001 \
  d005o000:005:0000 \
  d005p005:005:+005 \
  d006n001:006:-001 \
  d006n005:006:-005 \
  d006n006:006:-006 \
  d006p003:006:+003 \
  d006p006:006:+006 \
  d007o000:007:0000 \
  d007p003:007:+003 \
  d007p007:007:+007 \
  d012p007:012:+007 \
  d017p007:017:+007

RUNS_LATER :=

SHOW := 1
FORMAT := ps

.PHONY:: all single

all: 
	for f in ${RUNS}; do \
          ${MAKE} RUNARGS=$$f single; \
        done

# PSVIEW := kghostview
PSVIEW := evince

######################################################################
# Recursive make - caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS   := ${subst :, ,${RUNARGS}}
FUNCTAG  := ${word 1,${FIELDS}}
DEGREE   := ${word 2,${FIELDS}}
ORDER    := ${word 3,${FIELDS}}

FUNCDIR  := ${ALLFUNDIR}/${FUNCTAG}
OUTNAME  := ${FUNCDIR}/f
FUNCFILE := ${OUTNAME}.sfn

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

single: ${FUNCFILE} ${PSTARGET}

${FUNCFILE}: ${PROGDIR}/${MAKEPROG}
	mkdir -p ${FUNCDIR}
	${PROGDIR}/${MAKEPROG} \
          -degree ${DEGREE} -order ${ORDER} \
          -outName ${OUTNAME}

${PSTARGET}: \
            ${FUNCFILE} \
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
# ${RUNARGS}
######################################################################
