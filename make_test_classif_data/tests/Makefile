# Last edited on 2010-05-24 20:20:36 by stolfi

PROG := make_test_classif_data
PROGDIR := ..

VERIFY := YES

# PROBLEM:NA:NC:NS:SEED:NOISE:GRID

TESTS := \
  shells:02:04:000300:003173:0.00:Y \
  mballs:02:03:000300:003173:0.00:Y

TESTS_DONE := \
  saturn:02:02:000300:003173:0.00:N \
  saturn:02:02:003000:003173:0.00:Y \
  \
  petals:02:04:000300:003173:0.00:N \
  vessel:02:03:000300:003173:0.00:N \
  vessel:02:03:000300:003173:0.02:N \
  \
  mballs:02:03:000300:003173:0.00:N \
  shells:02:04:000300:003173:0.00:N \
  \
  shells:02:04:000300:003173:0.00:Y \
  mballs:02:03:000300:003173:0.00:Y \
  \
  mballs:02:03:010000:003173:0.00:N \
  mballs:02:03:010000:003173:0.00:Y
  
include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all:
	for t in ${TESTS}; do \
          ${MAKE} TEST=$$t single ; \
        done
        
clean::
	( cd out && rm -f *.dat *.ppm *.eps )
        
ifneq "/${TEST}" "/"        
########################################################################
# Make with ${TEST}

TEST_FIELDS := ${subst :, ,${TEST}}
PROBLEM := ${word 1, ${TEST_FIELDS}}
NA :=      ${word 2, ${TEST_FIELDS}}
NC :=      ${word 3, ${TEST_FIELDS}}
NS :=      ${word 4, ${TEST_FIELDS}}
SEED :=    ${word 5, ${TEST_FIELDS}}
NOISE :=   ${word 6, ${TEST_FIELDS}}
GRID :=    ${word 7, ${TEST_FIELDS}}

NAME :=  ${subst .,v,${subst :,-,${TEST}}}

single: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -problem ${PROBLEM} \
            -samples ${NS} \
            -seed ${SEED} \
            -prefix out/${NAME} \
            -attributes ${NA} 0 \
            -classes ${NC} \
            -noise ${NOISE} \
            -grid ${GRID} \
            -verify ${VERIFY} \
            -image 128 2 \
          > out/${NAME}.dat
	plot-dataset.sh ${NAME} ${NC}
	-display out/${NAME}.ppm
        
        
# end ${TEST}
########################################################################endif
endif
