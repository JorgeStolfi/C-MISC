# Last edited on 2017-02-26 03:36:40 by stolfilocal

PROG := SOMakeWaveBasis

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

MFREQS := \
  2:0002 \
  2:0020 \
  2:0200

# all: teste1
all: run-all

bas/wave/%.bas: ${PROGDIR}/SOMakeWaveBasis
	${PROGDIR}/SOMakeWaveBasis \
          -maxFreq $* \
          -dDim = 2 \
          -outName out/bas/wave/$*
