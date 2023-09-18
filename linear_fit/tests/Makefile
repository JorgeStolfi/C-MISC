# Last edited on 2023-04-26 22:06:29 by stolfi

PROG := linear_fit
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATA_TYPE}:{NX}:{WT_FLAG}:{UNIT_EXPL}:{UNIT_IMPL}:{VERBOSE}
# Were 
#    {NX} is the number of indep variable columns in the file.
#    {WT_FLAG} is 1 or 0 depending on whether there is a weight column in the file.
#    {UNIT_EXPL} is 1 if one of the {NX} variables is a constant "1".
#    {UNIT_IMPL} is "T" if the regression is to use an implicit "1" term, "F" otherwise. 
#    {VERBOSE} is 0 or 1.

EXACT_TESTS := \
  afffun:01:1:0:F:1 \
  afffun:01:1:0:T:1

RANDOM_TESTS := \
  rangen:04:0:0:F:1 \
  rangen:04:1:0:F:1 \
  rangen:05:0:1:F:1 \
  rangen:05:1:1:F:1 \
    \
  rangen:04:0:0:T:1 \
  rangen:04:1:0:T:1 \

FOC_TESTS := \
  fssqsq:03:0:0:F:1 \
  fssqsq:03:0:0:T:1
  
MORE_FOC_TESTS := \
  fshahi:02:0:0:F:1 \
  fshahi:02:0:0:T:1 \
   \
  fshahi:04:0:0:F:1 \
  fshahi:04:0:0:T:1 \
  fshahi:04:1:0:T:1 \
   \
  fshahi:03:0:0:F:1 \
  fshahi:03:0:0:T:1 \
   \
  fshahi:05:0:1:F:1 \
  fshahi:05:1:1:F:1 \
   \
  fshalo:05:1:1:F:1

ALL_TESTS := \
  fshasq:45:1:1:F:0 \
  fshaun:45:1:1:F:0 \
   \
  rangen:45:1:1:F:0 \
   \
  rangen:02:1:0:F:0 \
  rangen:04:0:0:F:0 \
  rangen:09:0:0:F:0 \
   \
  rangen:02:0:0:F:0 \
  rangen:04:1:0:F:0 \
  rangen:09:1:0:F:0 \
   \
  afffun:01:1:0:F:0 \
  bccbin:05:1:0:F:0 \
  bccfib:05:1:0:F:0 \
  bccpri:04:1:0:F:0 \
  bcfbin:05:1:0:F:0 \
  bcffib:05:1:0:F:0 \
  bcfpri:05:1:0:F:0 \
  bcgfib:04:1:0:F:0 \
  bcgpri:04:1:0:F:0 \
  brobro:01:1:0:F:0 \
  brobro:02:1:0:F:0 \
  expfun:02:1:0:F:0 \
  msnfun:04:1:0:F:0 \
  planou:03:0:1:F:0 \
  planou:03:1:1:F:0
 
# TESTS := ${ALL_TESTS}
# TESTS := ${RANDOM_TESTS}
# TESTS := ${FOC_TESTS}
TESTS := ${EXACT_TESTS}

.PHONY:: make-data do-all-tests do-single-test

all: make-data do-all-tests

do-all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-test ; \
        done

ifneq "/${TEST}" "/"
# ======================================================================
# Recursive section, client must define ${TEST}
  
TEST_ARGS := ${subst :, ,${TEST}}
DATA_TYPE :=  ${word 1,${TEST_ARGS}}
NX :=         ${word 2,${TEST_ARGS}}
WT_FLAG :=     ${word 3,${TEST_ARGS}}
UNIT_EXPL :=   ${word 4,${TEST_ARGS}}
UNIT_IMPL :=   ${word 5,${TEST_ARGS}}
VERBOSE :=    ${word 6,${TEST_ARGS}}

SHOW := SHOW

FULL_DATA_NAME := ${DATA_TYPE}_nx${NX}_wt${WT_FLAG}_ui${UNIT_EXPL}

IN_FILE := data/${FULL_DATA_NAME}.txt

OUT_PREFIX := out/${FULL_DATA_NAME}

OUT_FILE :=  ${OUT_PREFIX}.txt
FRM_FILE :=  ${OUT_PREFIX}.frm
PLOT_FILE := ${OUT_PREFIX}.png

TNAMES := \
  ${wordlist 1,${NX}, \
    AA AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU AV AW AX AY AZ \
    BA BB BC BD BE BF BG BH BI BJ BK BL BM BN BO BP BQ BR BS BT BU BV BW BX BY BZ \
    CA CB CC CD CE CF CG CH CI CJ CK CL CM CN CO CP CQ CR CS CT CU CV CW CX CY CZ \
    DA DB DC DD DE DF DG DH DI DJ DK DL DM DN DO DP DQ DR DS DT DU DV DW DX DY DZ \
  }
  
TITLE := "Test of {linear_fit} - data ${DATA_TYPE}_nx${NX} - weight ${WT_FLAG} - unit expl ${UNIT_EXPL} impl ${UNIT_IMPL}"
TITLE_DATA := "Given"
TITLE_FIT := "Fitted"
        
ifeq "/${DATA_TYPE}" "/rangen"
  # Randomly generated dataset:
  MAKER := make_regr_data_file.gawk
  ${IN_FILE}: ${MAKER}
	  ./${MAKER} -v fname="${FULL_DATA_NAME}" > ${IN_FILE}
endif

do-single-test: ${PROGDIR}/${PROG} data/${FULL_DATA_NAME}.txt 
	rm -f ${OUT_PREFIX}.{txt,png}
	${PROGDIR}/${PROG} \
          -terms ${NX} \
          -termNames ${TNAMES} \
          -weighted ${WT_FLAG} \
          -verbose ${VERBOSE} \
          -writeFormula ${FRM_FILE} \
          -unitTerm ${UNIT_IMPL} \
          -format '%+14.8f' \
          < ${IN_FILE} \
          > ${OUT_FILE}
	plot_regression_result.sh ${SHOW} ${OUT_PREFIX} ${TITLE} ${TITLE_DATA} ${TITLE_FIT}
        
# ======================================================================
endif
