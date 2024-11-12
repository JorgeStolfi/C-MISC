# Last edited on 2008-06-10 20:42:39 by stolfi

TESTPROG := SPHessTest

PROGS := ${TESTPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

# FNAME:GNAME:TRINAME:SMPORDER:EVDOT:VGDOT:SLDOT

# The integral of {x^n} over the sphere is 
# {4*PI/(n+1)} if {n} is even, 0 if {n} is odd.

RUNS := \
  proc/sepx:proc/linx:icosa/raw:40:1.39626340159546:0.00000000000000:0.00000000000000

RUNS_DONE := \
  harm/d001n001:harm/d001n001:icosa/raw:40:1.00000000000000:0.00000000000000:0.00000000000000 \
  harm/d001o000:harm/d001o000:icosa/raw:40:1.00000000000000:0.00000000000000:0.00000000000000 \
   \
  harm/d000o000:harm/d001o000:icosa/raw:40:0.00000000000000:0.00000000000000:0.00000000000000 \
  harm/d001n001:harm/d001o000:icosa/raw:40:0.00000000000000:0.00000000000000:0.00000000000000 \
  harm/d006n006:harm/d001o000:icosa/raw:40:0.00000000000000:0.00000000000000:0.00000000000000 \
   \
  proc/linx:proc/linx:icosa/raw:40:4.18879020478639:0.00000000000000:0.00000000000000 \
  proc/sqrx:proc/linx:icosa/raw:40:0.00000000000000:0.00000000000000:0.00000000000000 \
   \
  proc/hrm2:proc/unit:icosa/raw:40:0.00000000000000:0.00000000000000:0.00000000000000
  

.PHONY:: single

all: 
	for r in ${RUNS}; do \
          ${MAKE} RUNARGS=$$r single; \
        done

# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS   := ${subst :, ,${RUNARGS}}
FFUNC    := ${word 1,${FIELDS}}
GFUNC    := ${word 2,${FIELDS}}
TRIANG   := ${word 3,${FIELDS}}
SMPORDER := ${word 4,${FIELDS}}
EVDOT    := ${word 5,${FIELDS}}
VGDOT    := ${word 6,${FIELDS}}
SLDOT    := ${word 7,${FIELDS}}

FNAME := ${DATADIR}/fn/${FFUNC}/f
GNAME := ${DATADIR}/fn/${GFUNC}/f
TRINAME := ${DATADIR}/sp/pw/${TRIANG}/t

TESTPROG := SPDotTest

single: \
            ${FNAME}.sfn ${GNAME}.sfn ${TRINAME}.tri \
            ${PROGDIR}/${TESTPROG} \
            Makefile
	${PROGDIR}/${TESTPROG} \
          -funcNames ${FNAME} ${GNAME} \
          -triName ${TRINAME} \
          -smpOrder ${SMPORDER} \
          -exact ${EVDOT} ${VGDOT} ${SLDOT}

endif
# ${RUNARGS}
######################################################################
