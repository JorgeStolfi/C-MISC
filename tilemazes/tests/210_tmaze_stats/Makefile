# Last edited on 2023-02-04 07:02:34 by stolfi

PROGDIR := ../..
PROG := tmaze_stats

# ${FAMILY}-${NX}-${NY}-${SEED}-${TOP_TAG}-${TRIALS}

SMALL_TESTS := \
  brasilia-000007-000005-0000002003-o-000002 \
  brasilia-000007-000005-0000000001-t-000002 \
   \
  blipblop-000007-000005-0000002003-o-000002 \
  blipblop-000007-000005-0000000001-t-000002 \

MEDIUM_TESTS := \
  brasilia-000200-000200-0000002003-o-000010 \
  brasilia-000200-000200-0000000001-t-000010 \
   \
  blipblop-000200-000200-0000002003-o-000010 \
  blipblop-000200-000200-0000000001-t-000010 \

BIG_TESTS := \
  brasilia-001000-001000-0000004615-t-000100 \
  blipblop-001000-001000-0000004615-t-000100 \
   \
  brasilia-002000-002000-0000004651-t-000100 \
  blipblop-002000-002000-0000004651-t-000100 \
   \
  brasilia-003000-003000-0000004643-t-000100 \
  blipblop-003000-003000-0000004643-t-000100
  
TESTS := ${SMALL_TESTS} ${MEDIUM_TESTS}
# TESTS := ${BIG_TESTS}

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
TRIALS :=    ${word 6,${TESTFIELDS}}

TOPOLOGY := ${subst /o,open,${subst /t,torus,/${TOP_TAG}}}

OUTNAME := out/${FAMILY}/${NX}-${NY}-${SEED}-${TOP_TAG}-${TRIALS}

single-test: ${PROGDIR}/${PROG}
	-rm -f ${OUTNAME}.*
	-${PROGDIR}/${PROG} \
          -family ${FAMILY} \
          -size ${NX} ${NY} \
          -seed ${SEED} \
          -topology ${TOPOLOGY} \
          -trials ${TRIALS}
	-if [[ -s ${OUTNAME}-small.act ]] ; then plot-comp-sizes.sh ${OUTNAME}-small.act ; fi
	-if [[ -s ${OUTNAME}-large.act ]] ; then plot-comp-sizes.sh ${OUTNAME}-large.act ; fi

endif
# end ${TESTARGS}
######################################################################
