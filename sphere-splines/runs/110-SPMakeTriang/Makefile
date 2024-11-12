# Last edited on 2023-02-04 07:04:05 by stolfi

PROGS := \
  SPMakeTriang \
  SPPlotTriang
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

RUNS := \
  octa:raw \
  icosa:raw
  
NONRUNS :=  \
  tetra:raw

.PHONY:: \
  run clean-single run-single show-single

all: clean run

clean::
	for d in ${RUNS} ; do \
          ${MAKE} RUNARGS=$$d clean-single ; \
        done
  
run: 
	for d in ${RUNS} ; do \
          ${MAKE} RUNARGS=$$d run-single show-single ; \
        done

# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"
FIELDS  := ${subst :, ,${RUNARGS}}
TRI    := ${word 1,${FIELDS}}
REFRES := ${word 2,${FIELDS}}

TRIDIR := ${DATADIR}/sp/pw/${TRI}/${REFRES}
TRINAME := ${TRIDIR}/t
PSNAME := ${TRINAME}-py

clean-single:
	/bin/rm -f ${TRINAME}.tri

run-single: ${TRINAME}.tri

show-single: ${PSNAME}.ps

${TRIDIR}: 
	mkdir -p ${TRIDIR}
  
# PSVIEW := kghostview
PSVIEW := evince

# PAPEROPTIONS := -figSize 230 -paperSize a3
PAPEROPTIONS := -figSize 150 -paperSize letter

${PSNAME}.ps: ${PROGDIR}/SPPlotTriang Makefile \
            ${TRIDIR} ${TRINAME}.tri
	${PROGDIR}/SPPlotTriang \
          -triName ${TRINAME} \
          -outName ${PSNAME} \
          -obs 0 4 2 1 -radius 1.05 \
          -showBack ${PAPEROPTIONS} \
          -caption "${TRINAME}.tri ${RUNARGS}"
	${PSVIEW} ${PSNAME}.ps

# ======================================================================
# Rules for individual triangulations

${DATADIR}/sp/pw/tetra/raw/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -tetrahedron

${DATADIR}/sp/pw/octa/raw/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	@echo "LD_LIBRRAY_PATH=$${LD_LIBRRAY_PATH}"
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -octahedron

${DATADIR}/sp/pw/icosa/raw/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -icosahedron

${DATADIR}/sp/pw/rand/r0100/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -nSites 100

${DATADIR}/sp/pw/rand/r0200/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -nSites 200

${DATADIR}/sp/pw/rand/r0400/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -nSites 400

${DATADIR}/sp/pw/rand/r0800/t.tri: ${PROGDIR}/SPMakeTriang Makefile \
            ${TRIDIR}
	${PROGDIR}/SPMakeTriang \
          -outName ${TRINAME} \
          -nSites 800

endif
# End of ${RUNARGS} part.
# ======================================================================
