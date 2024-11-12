# Last edited on 2023-02-04 07:03:59 by stolfi

PROGS := \
  SPMake3Triang \
  SPPlotTriang
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

RUNS := \
  icosa:raw:r3/01 \
  icosa:r3/01:r3/02 \
  icosa:r3/02:r3/03 \
  icosa:r3/03:r3/04 \
  icosa:r3/04:r3/05
  
NONRUNS :=

.PHONY:: run clean-single run-single

all: run

clean::
	for tri in ${RUNS} ; do \
          ${MAKE}  RUNARGS="$$tri" clean-single ; \
        done
  
run: 
	for tri in ${RUNS} ; do \
          ${MAKE}  RUNARGS="$$tri" run-single ; \
        done
        
# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS   := ${subst :, ,${RUNARGS}}
TRI      := ${word 1,${FIELDS}}
TAGIN    := ${word 2,${FIELDS}}
TAGOT    := ${word 3,${FIELDS}}

INDIR := ${DATADIR}/sp/pw/${TRI}/${TAGIN}
INNAME := ${INDIR}/t

OTDIR := ${DATADIR}/sp/pw/${TRI}/${TAGOT}
OTNAME := ${OTDIR}/t
PSNAME := ${OTNAME}-py

clean-single: 
	-/bin/rm -f ${OTNAME}.tri ${OTNAME}.log  ${PSNAME}.ps
  
run-single: clean-single ${OTNAME}.tri ${PSNAME}.ps

${OTDIR}:
	mkdir -p ${OTDIR}

${OTNAME}.tri: ${OTDIR} ${INNAME}.tri \
            ${PROGDIR}/SPMake3Triang Makefile
	${PROGDIR}/SPMake3Triang \
            -triName ${INNAME} \
            -outName ${OTNAME} \
          2>&1 | tee ${OTNAME}.log

#  -verbose -plot \
# \
#

# PSVIEW := kghostview
PSVIEW := evince

# PAPEROPTIONS := -figSize 230 -paperSize a3
PAPEROPTIONS := -figSize 150 -paperSize letter

${PSNAME}.ps: ${PROGDIR}/SPPlotTriang Makefile \
            ${OTNAME}.tri
	${PROGDIR}/SPPlotTriang \
          -triName ${OTNAME} \
          -outName ${PSNAME} \
          -obs 0 4 2 1 -radius 1.05 \
          -showBack ${PAPEROPTIONS} \
          -caption "${OTNAME}.tri ${RUNARGS}"
	${PSVIEW} ${PSNAME}.ps

endif
# End of ${RUNARGS} part.
# ======================================================================
