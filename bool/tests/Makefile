# Last edited on 2023-02-08 21:52:32 by stolfi

PROG := bool

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

SHELL =	/bin/bash

.PHONY:: all-tests

TESTS := \
  1-2:0010 \
  2-1:0100 \
  1+2:0111 \
  2+1:0111 \
  1.2:0001 \
  2.1:0001 \
  1~2:0110 \
  2~1:0110
  
INPUTS := A.txt B.txt

SUBSETS := R-0001.txt R-0010.txt R-0100.txt R-0110.txt R-0111.txt

all: all-tests

all-tests: 
	for tt in ${TESTS} ; do \
          for sgn in s u ; do \
            ${MAKE} TEST="$$tt:$$sgn" single ; \
          done ; \
        done
        
# ######################################################################
# Recursive section - needs ${TEST}
ifneq "/${TEST}" "/"

TARGS := ${subst :, , ${TEST}}
OP :=  ${word 1,${TARGS}}
RES := ${word 2,${TARGS}}
SGN := ${word 3,${TARGS}}

IN_FILES := data/A-${SGN}.txt data/B-${SGN}.txt
REF_FILE := data/R-${RES}-${SGN}.txt

single: ${IN_FILES} ${REF_FILE}
	bool -${SGN} "${OP}" ${IN_FILES} > .aaa
	diff .aaa ${REF_FILE} > .bbb
	@if [[ -s .bbb ]] ; then echo "** FAILED"; cat .bbb ; else echo "~~~ OK! ~~~"; fi

endif
# End ${TEST} section
# ######################################################################
clean::

