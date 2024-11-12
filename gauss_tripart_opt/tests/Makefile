# Last edited on 2020-01-15 19:31:13 by jstolfi

PROG := gauss_tripart_opt
PROGDIR := ..

.PHONY:: single

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

HALFF := 0.5
SQRTH := 0.70710678118654752440

TESTS := \
  sqrth_1_F:1:F:0,0,0,${SQRTH},${SQRTH}/1,0.5000,1.5000,${SQRTH},${SQRTH} \
  halff_3_F:3:F:0,0,0,${HALFF},${HALFF}/1,0.4500,0.5500,${HALFF},${HALFF}/2,0.9500,1.1000,${HALFF},${HALFF},${HALFF}/3,1.4500,1.5500,${HALFF},${HALFF} \

MORE_TESTS := \
  sqrth_1_T:1:T:0,0,0,${SQRTH},${SQRTH}/1,0.5000,1.5000,${SQRTH},${SQRTH} \
  halff_2_F:2:F:0,0,0,${HALFF},${HALFF}/1,0.3000,0.7000,${HALFF},${HALFF}/2,0.8000,1.5000,${HALFF},${HALFF} \
  halff_1_F:1:F:0,0,0,${HALFF},${HALFF}/1,0.3000,1.5000,${HALFF},${HALFF} \
  halff_1_T:1:T:0,0,0,${HALFF},${HALFF}/1,0.3000,1.5000,${HALFF},${HALFF} \
  halff_2_F:2:F:0,0,0,${HALFF},${HALFF}/1,0.3000,0.7000,${HALFF},${HALFF}/2,0.8000,1.5000,${HALFF},${HALFF} \
  
all: ${PROGDIR}/${PROG} do_test.sh Makefile
	for tst in ${TESTS}; do \
          ${MAKE} ARGS="$${tst}" single ; \
        done
  
ifneq "/${ARGS}" "/"
ARGWORDS = ${subst :, ,${ARGS}}
TAG      := ${word 1,${ARGWORDS}}
NTERMS   := ${word 2,${ARGWORDS}}
NOMIDDLE := ${word 3,${ARGWORDS}}
RANGES   := ${word 4,${ARGWORDS}}

single: ${PROGDIR}/${PROG} do_test.sh 
	-rm -f out/funcs_${TAG}.*
	do_test.sh ${TAG} ${NTERMS} ${NOMIDDLE} "${RANGES}"
	${PROGDIR}/plot_gto_funcs.sh ${TAG} ${NTERMS}
 
endif
