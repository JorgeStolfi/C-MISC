# Last edited on 2023-02-21 10:42:45 by stolfi

PROG := ioplot
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: check

check: \
  run-test-1 \
  run-test-2 \
  run-test-3
  
PSVIEW := evince
PNGVIEW := display
  
# ====================================================================
# Tests:

OUTDIR := out

run-test-1: ${PROGDIR}/${PROG} ${OUTDIR}
	${PROGDIR}/${PROG} \
            -figSize 180 150 \
            -plotWidth 80 \
            -debug \
          < testfile.dat \
          > ${OUTDIR}/testfile-1.eps
	${PSVIEW} ${OUTDIR}/testfile-1.eps
	convert -density 150 -background white ${OUTDIR}/testfile-1.eps -flatten ${OUTDIR}/testfile-1.png
	${PNGVIEW} ${OUTDIR}/testfile-1.png

run-test-2: ${PROGDIR}/${PROG} ${OUTDIR}
	${PROGDIR}/${PROG} \
            -figSize 180 150 \
            -plotWidth 80 \
            -simplify \
            -debug \
          < testfile.dat \
          > ${OUTDIR}/testfile-2.eps
	${PSVIEW} ${OUTDIR}/testfile-2.eps
	convert -density 150 -background white ${OUTDIR}/testfile-2.eps -flatten ${OUTDIR}/testfile-2.png 
	${PNGVIEW} ${OUTDIR}/testfile-2.png

run-test-3: ${PROGDIR}/${PROG} ${OUTDIR}
	${PROGDIR}/${PROG} \
            -figSize 180 150 \
            -plotWidth 80 \
            -showNodes \
            -font "TimesRoman" 12 \
            -fillColor 0.400 0.600 1.000 \
            -drawColor 0.000 0.000 0.500 \
            -textColor 0.000 0.000 0.500 \
            -depth all 1 \
            -debug \
          < testfile.dat \
          > ${OUTDIR}/testfile-3.eps
	${PSVIEW} ${OUTDIR}/testfile-3.eps
	convert -density 150 -background white ${OUTDIR}/testfile-3.eps -flatten ${OUTDIR}/testfile-3.png
	${PNGVIEW} ${OUTDIR}/testfile-3.png 
  
${OUTDIR}:
	mkdir -p ${OUTDIR}

clean::
	rm -f ${OUTDIR}/*.eps
