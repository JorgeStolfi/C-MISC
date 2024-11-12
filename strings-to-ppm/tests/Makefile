# Last edited on 2012-07-22 14:34:56 by stolfilocal

PROG := ../strings-to-ppm

TESTNAME := test1

SHOWPPM := display

PPMFILE := ${TESTNAME}.ppm
STRINGFILE := ${TESTNAME}.str

all: ${PPMFILE}

${PPMFILE}: ${STRINGFILE} ${PROG}
	${PROG} -numBits 16 < ${STRINGFILE} > ${PPMFILE}
	${SHOWPPM} ${PPMFILE}

clean::

actions:
	${MAKE} ${ACTIONS}

