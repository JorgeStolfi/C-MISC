# Makefile for testing the wavefront simulator
# Last edited on 2023-02-04 07:04:17 by stolfi

PROGS := \
  quadwfc \
  plotwfc \
  povwfc

PROGDIR := ..

.PHONY:: \
  run-tests run-test-quadwfc run-test-plotwfc run-test-povwfc 

all: run-tests

# run-tests: run-test-quadwfc run-test-plotwfc
run-tests: run-test-quadwfc run-test-povwfc
# run-tests: run-test-plotwfc

clean::

actions:
	${MAKE} ${ACTIONS}

# ====================================================================
# Running tests

TNAME := testfile
OUTPREFIX := out/${TNAME}

MAXITER := 30

# FRAME := ini
# FRAME := 000002
# FRAME := 000003
# FRAME := 000004
FRAME := 000010

# Simulator:

TRIFILE := ${OUTPREFIX}-${FRAME}.tri

run-test-quadwfc: ${TRIFILE}

${TRIFILE}: ${PROGDIR}/quadwfc ${TNAME}.dat Makefile
	${PROGDIR}/quadwfc \
            -outName ${OUTPREFIX} \
            -maxIter ${MAXITER} \
            -timeStep 5.0 \
            -edgeLength 200 \
            -signature P0rP \
          < ${TNAME}.dat

# Postscript plotting:

EPSFILE := ${OUTPREFIX}-${FRAME}.eps

PSVIEW := evince

run-test-plotwfc: ${EPSFILE}

${EPSFILE}: ${PROGDIR}/plotwfc ${TRIFILE} Makefile
	${PROGDIR}/plotwfc \
            -triName ${OUTPREFIX}-${FRAME} \
            -geoName ${TNAME} \
            -obs 0 1 2 3 \
            -ctr 1 5000 5000 3000 \
            -up 0 0 0 1 \
            -radius 5000 -eps \
            -outName ${OUTPREFIX}-${FRAME}
	${PSVIEW} ${EPSFILE}

# POV-Ray plotting:

INCFILE := ${OUTPREFIX}-${FRAME}.inc
PNGFILE := ${OUTPREFIX}-${FRAME}.png

run-test-povwfc: ${PNGFILE}

${INCFILE}: ${PROGDIR}/povwfc ${TRIFILE} Makefile
	${PROGDIR}/povwfc \
            -triName ${OUTPREFIX}-${FRAME} \
            -geoName ${TNAME} \
            -outName ${OUTPREFIX}-${FRAME} \
            -lineWidth 20
            
OTHERINCS := \
  pov-main.pov \
  pov-lights-camera.inc \
  pov-textures.inc

${PNGFILE}: pov-render ${INCFILE} ${OTHERINCS} Makefile
	pov-render -show ${INCFILE}
