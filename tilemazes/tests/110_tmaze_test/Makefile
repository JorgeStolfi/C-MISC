# Last edited on 2023-02-04 07:00:02 by stolfi

PROGDIR := ../..
PROG := tmaze_test

# ${FAMILY}-${NX}-${NY}-${SEED}-${TOP_TAG}-${CELL_SIZE}

TESTS := \
  brasilia-000003-000002-0000000001-o-10.000 \
  brasilia-000039-000054-0000000001-t-04.233 \
   \
  blipblop-000003-000002-0000000001-o-10.000 \
  blipblop-000039-000054-0000000001-t-04.233

.PHONY:: all-tests single-test

all: all-tests

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all-tests:
	for tt in ${TESTS}; do \
          ${MAKE} TESTARGS="$$tt" single-test ; \
        done

# PSVIEW := gv
PSVIEW := evince

######################################################################
# For recursive make - caller must define ${TESTARGS}

ifneq "/${TESTARGS}" "/"

TESTFIELDS := ${subst -, ,${TESTARGS}}
FAMILY :=    ${word 1,${TESTFIELDS}}
NX :=        ${word 2,${TESTFIELDS}}
NY :=        ${word 3,${TESTFIELDS}}
SEED :=      ${word 4,${TESTFIELDS}}
TOP_TAG :=   ${word 5,${TESTFIELDS}}
CELL_SIZE := ${word 6,${TESTFIELDS}}

TOPOLOGY := ${subst /o,open,${subst /t,torus,/${TOP_TAG}}}

OUTNAME := out/${FAMILY}/${NX}-${NY}-${SEED}-${TOP_TAG}

single-test: ${PROGDIR}/${PROG}
	-rm -f ${OUTNAME}.*
	-${PROGDIR}/${PROG} \
          -family ${FAMILY} \
          -size ${NX} ${NY} \
          -seed ${SEED} \
          -topology ${TOPOLOGY} \
          -print \
          -graph \
          -comps \
          -plot ${CELL_SIZE} -cells -grid -style curved
	-if [[ -s ${OUTNAME}.prt ]] ; then cat ${OUTNAME}.prt ; fi
	-if [[ -s ${OUTNAME}.cts ]] ; then cat ${OUTNAME}.cts ; fi
	-if [[ -s ${OUTNAME}.eps ]] ; then ${PSVIEW} ${OUTNAME}.eps ; fi

endif
# end ${TESTARGS}
######################################################################
