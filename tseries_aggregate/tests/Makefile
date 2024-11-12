# Last edited on 2014-01-14 20:17:46 by stolfilocal

PROG := tseries_aggregate
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Format: {DATASET}:{NG}:{GROUPING}

TESTS := \
  exp:2:none \
  bro:1:none \
  bro:2:none
 
ALL_TESTS := \
  seq:0:none \
  seq:1:none \
   \
  seq:5:none \
  seq:5:binary \
  seq:5:fibonacci \
   \
  gox:4:none \
  gox:4:fibonacci \
   \
  bcc:4:none \
   \
  bcc:5:binary \
  bcc:5:fibonacci \
   \
  bro:1:none \
  bro:2:none \
   \
  lin:1:none \
  sin:2:none \
  msn:4:none \
  exp:2:none
  

.PHONY:: do-all-tests do-single-test

all: do-all-tests

do-all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt do-single-test ; \
        done

ifneq "/${TEST}" "/"
# ======================================================================
# Recursive section, client must define ${TEST}
  
TEST_ARGS := ${subst :, ,${TEST}}
DATASET :=  ${word 1,${TEST_ARGS}}
NG :=       ${word 2,${TEST_ARGS}}
GROUPING := ${word 3,${TEST_ARGS}}

SHOW := SHOW

IN_FILE := data/${DATASET}.txt

OUT_PREFIX := out/${DATASET}_n${NG}_${GROUPING}

OUT_FILE := ${OUT_PREFIX}.txt

PLOT_FILE := ${OUT_PREFIX}.png

do-single-test: ${PROGDIR}/${PROG} ${IN_FILE}
	rm -f ${OUT_PREFIX}.{txt,png}
	${PROGDIR}/${PROG} \
          -terms ${NG} \
          -group ${GROUPING} \
          -constant 1 \
          < ${IN_FILE} \
          > ${OUT_FILE}
        
# ----------------------------------------------------------------------
# Obtaining the test data:

MKREC := make_recursive_test_data.gawk
MKSIN := make_multisin_test_data.gawk

data/lin.txt: ${MKREC}
	 ${MKREC} -v NX=3 -v C="1,1" -v X="0" -v D="0" > data/lin.txt

data/exp.txt: ${MKREC}
	 ${MKREC} -v NX=3000 -v C="0,1.001" -v X="1" -v D="0" > data/exp.txt

data/bro.txt: ${MKREC}
	 ${MKREC} -v NX=3000 -v C="0,1.000" -v X="0" -v D="1.0" > data/bro.txt

data/sin.txt: ${MKSIN}
	 ${MKSIN} -v NX=3000 -v A="4" -v F="9" -v P="0" > data/sin.txt

data/msn.txt: ${MKSIN}
	 ${MKSIN} -v NX=3000 -v A="3,2,1,1" -v F="4,7,9,5" -v P="0.0,0.5,0.3,0.7" > data/msn.txt

# ======================================================================
endif
