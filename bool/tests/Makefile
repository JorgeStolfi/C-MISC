# Last edited on 2008-01-15 00:40:04 by stolfi

PROG := bool

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

SHELL =	/bin/sh

.PHONY:: all-tests

TESTS := \
  test-d12 \
  test-d21 \
  test-u12 \
  test-i12 \
  test-s12
  
INPUTS := A.txt B.txt

SUBSETS := R-0001.txt R-0010.txt R-0100.txt R-0110.txt R-0111.txt

all: all-tests

all-tests: ${TESTS}

test-d12: ${INPUTS} ${SUBSETS}
	bool '1-2' A.txt B.txt > .aaa
	diff .aaa R-0010.txt

test-d21: ${INPUTS} ${SUBSETS}
	bool '2-1' A.txt B.txt > .aaa
	diff .aaa R-0100.txt

test-u12: ${INPUTS} ${SUBSETS}
	bool '2+1' A.txt B.txt > .aaa
	diff .aaa R-0111.txt

test-i12: ${INPUTS} ${SUBSETS}
	bool '2.1' A.txt B.txt > .aaa
	diff .aaa R-0001.txt

test-s12: ${INPUTS} ${SUBSETS}
	bool '2#1' A.txt B.txt > .aaa
	diff .aaa R-0110.txt

clean::
