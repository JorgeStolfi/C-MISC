# Last edited on 2012-07-22 14:34:31 by stolfilocal
  
PROG := sort-distr
  
PROGDIR := ..

RUNS := \
  cr-eu-pr/circle \
  pt-eu-pr/path \
  mw-eu-pr/weight \
  mw-eu-ds/weight \
  dj-eu-pr/geom \
  dj-gm-pr/geom
  
.SUFFIXES:
.PHONY:: single

NONRUNARGS :=

SHOWPPM := display -filter Box -resize '1000%!'
  
all: 
	for t in ${RUNS}; do \
          ${MAKE} RUNARGS=$$t single ; \
        done

clean::

actions:
	${MAKE} ${ACTIONS}

######################################################################
# Recursive make - caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

ARGS := ${subst /, ,${RUNARGS}}
TESTNAME  := ${word 1,${ARGS}}
INPNAME   := ${word 2,${ARGS}}

OPTFILE := tst-${TESTNAME}.args
INPFILE := tst-${INPNAME}.inp
TITFILE := tst-${TESTNAME}.tit
PPMFILE := tst-${TESTNAME}.ppm

single: ${OPTFILE} ${INPFILE} ${TITFILE} ${PROGDIR}/${PROG}
	@echo '====================================='
	@echo "input = ${INPNAME}  options = `cat ${TITFILE}`"
	${PROGDIR}/${PROG} -p ${PPMFILE} `cat ${OPTFILE}` < ${INPFILE}
	${SHOWPPM} ${PPMFILE}

endif
# end make ${RUNARGS}
######################################################################
