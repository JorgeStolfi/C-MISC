# Last edited on 2023-02-04 07:03:34 by stolfi

MAKEPROG := SPMakeTimeOnlyProcFunction
PLOTPROG := SPPlotTimeOnlyFunction

PROGS := ${MAKEPROG} ${PLOTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data
ALLFUNDIR := ${DATADIR}/tofn/proc

FUNS := \
  unit \
  lint \
  htsq \
  cost \
  exnt \
  amor
  
FUNS_LATER :=

SHOW := 0
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
FUNFILE := ${OUTNAME}.tof

single: ${FUNFILE}

${FUNFILE}: ${PROGDIR}/${MAKEPROG}
	mkdir -p ${FUNDIR}
	${PROGDIR}/${MAKEPROG} \
          -function ${FUN} \
          -outName ${OUTNAME}

endif
# ${FUN}
######################################################################
