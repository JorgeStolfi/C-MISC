# Last edited on 2017-02-26 03:35:30 by stolfilocal

PROG := SOMakeProcFunction

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

FUNS := \
  unit \
  linx \
  liny \
  linz \
  lind \
  sqrx \
  para \
  hump \
  sinx \
  sxcy \
  sxsy \
  gaus \
  dgdx \
  mhat \
  cbsm

FUNFILES := ${addsuffix .fun,${addprefix fun/,${FUNS}}}

all: ${FUNFILES}

fun/%.fun: ${PROGDIR}/SOMakeProcFunction
	${PROGDIR}/SOMakeProcFunction \
          -funcName $* \
          -pDim 2 -fDim 1 \
          -outName out/fun/$*
