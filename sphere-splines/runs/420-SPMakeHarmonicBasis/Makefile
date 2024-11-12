# Last edited on 2008-06-10 20:42:57 by stolfi

PROG := SPMakeHarmonicBasis

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

# CHECKOPT := -checkOrtho
CHECKOPT := 

# DEGREE:SMPORDER

RUNS :=

RUNS_TO_DO :=
  
RUNS_DONE := \
  004:40 \
  005:40 \
  006:40 \
  016:40 \
  017:40 \
  018:40 \
  012:40

.PHONY:: \
  run \
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
DEGREE    := ${word 1,${FIELDS}}
SMPORDER  := ${word 2,${FIELDS}}

BASDIR := ${DATADIR}/sp/gl/harm/d${DEGREE}
BASNAME := ${BASDIR}/b

run-single: ${BASNAME}.bas

${BASNAME}.bas: ${PROGDIR}/SPMakeHarmonicBasis
	mkdir -p ${BASDIR}
	${PROGDIR}/SPMakeHarmonicBasis \
          -degree ${DEGREE} \
          ${CHECKOPT} -verbose \
          -smpOrder ${SMPORDER} \
          -outName ${BASNAME}

endif
# end ${RUNARGS}
######################################################################
