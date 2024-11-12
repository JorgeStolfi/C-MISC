# Last edited on 2023-02-04 07:02:50 by stolfi

PROG := SPUniSolve

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 40

RUNS := \
  Fsepx:+0.500:sepx:pw/icosa/raw/h-d7c1 \
  Fsepx:+0.500:sepx:pw/icosa/raw/nh-d6c1 \
  \
  Foctx:+0.500:octx:pw/icosa/raw/h-d7c1 \
  Foctx:+0.500:octx:pw/icosa/raw/nh-d6c1 \
  \
  Fmcos:+0.500:cosx:pw/icosa/raw/h-d7c1 \
  Fmcos:+0.500:cosx:pw/icosa/raw/nh-d6c1 \
   \
  Fsinx:+0.500:sinx:pw/icosa/raw/h-d7c1 \
  Fsinx:+0.500:sinx:pw/icosa/raw/nh-d6c1 \
   \
  Fexpx:+0.500:expx:pw/icosa/raw/h-d7c1 \
  Fexpx:+0.500:expx:pw/icosa/raw/nh-d6c1 \
   \
  Fmexp:+0.500:expx:pw/icosa/raw/h-d7c1 \
  Fmexp:+0.500:expx:pw/icosa/raw/nh-d6c1 \

RUNS_DONE := \
  Fsqrx:+0.500:sqrx:pw/icosa/raw/h-d7c1 \
  Fsqrx:+0.500:sqrx:pw/icosa/raw/nh-d6c1 \
   \
  Fcubx:+0.500:cubx:pw/icosa/raw/h-d7c1 \
  Fcubx:+0.500:cubx:pw/icosa/raw/nh-d6c1 \

RUNS_MONDAY :=

RAW_ICOSA_RUNS_H7_NH6 := \
  Fmoct:+0.500:octx:pw/icosa/raw/h-d7c1 \
  Fmoct:+0.500:octx:pw/icosa/raw/nh-d6c1 \
   \
  Fexpx:+0.500:expx:pw/icosa/raw/h-d7c1 \
  Fexpx:+0.500:expx:pw/icosa/raw/nh-d6c1 \
  \
  Fmexp:+0.500:expx:pw/icosa/raw/h-d7c1 \
  Fmexp:+0.500:expx:pw/icosa/raw/nh-d6c1 \
  \
  Fcosx:+0.500:cosx:pw/icosa/raw/h-d7c1 \
  Fcosx:+0.500:cosx:pw/icosa/raw/nh-d6c1 \
  \
  Fsinx:+0.500:sinx:pw/icosa/raw/h-d7c1 \
  Fsinx:+0.500:sinx:pw/icosa/raw/nh-d6c1 \
  \
  Foctx:+0.500:octx:pw/icosa/raw/h-d7c1 \
  Foctx:+0.500:octx:pw/icosa/raw/nh-d6c1 \
  \
  Fqrtx:+0.500:qrtx:pw/icosa/raw/h-d7c1 \
  Fqrtx:+0.500:qrtx:pw/icosa/raw/nh-d6c1 \
  \
  Fmsep:+0.500:sepx:pw/icosa/raw/h-d7c1 \
  Fmsep:+0.500:sepx:pw/icosa/raw/nh-d6c1 \
   \
  Fmoct:+0.500:octx:pw/icosa/raw/h-d7c1 \
  Fmoct:+0.500:octx:pw/icosa/raw/nh-d6c1 \
   \
  Fcubx:+0.500:cubx:pw/icosa/raw/h-d7c1 \
  Fcubx:+0.500:cubx:pw/icosa/raw/nh-d6c1 \
  \
  Fsqrx:+0.500:sqrx:pw/icosa/raw/h-d7c1 \
  Fsqrx:+0.500:sqrx:pw/icosa/raw/nh-d6c1

OTHER_RUNS := \
  Fmexp:+0.500:expx:pw/icosa/r3/01/nh-d6c1 \
  Fmexp:+0.500:expx:pw/icosa/r3/01/h-d7c1 \
  \
  Fexpx:+0.500:expx:pw/icosa/raw/h-d5c1 \
  Fexpx:+0.500:expx:pw/icosa/raw/h-d6c1 \
  Fexpx:+0.500:expx:pw/icosa/raw/h-d7c1 \
  Fexpx:+0.500:expx:pw/icosa/raw/nh-d6c1 \
   \
  Fexpx:+0.500:expx:gl/bez/nh-d3 \
  Fexpx:+0.500:expx:gl/bez/nh-d4 \
  Fexpx:+0.500:expx:gl/bez/h-d4 \
   \
  Fqrtx:+0.500:qrtx:gl/bez/nh-d4 \
  Fmexp:+0.500:expx:gl/bez/nh-d4
  
BUGRUNS := \
  Fcosx:+0.500:cosx:gl/bez/nh-d4 \
  Fsinx:+0.500:sinx:gl/bez/nh-d4

.PHONY:: run clean-single run-single \
  make-approx show-ps

all: clean run
# all: run

clean::
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" clean-single ; \
        done

run: 
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" run-single ; \
        done
  
# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS  := ${subst :, ,${RUNARGS}}
RHSFUNC := ${word 1,${FIELDS}}
COEFF   := ${word 2,${FIELDS}}
SOLNAME := ${word 3,${FIELDS}}
BASIS   := ${word 4,${FIELDS}}

DATADIR := ../../data

BASDIR  := ${DATADIR}/sp/${BASIS}
BASNAME := ${BASDIR}/b

EVMATNAME := ${BASNAME}-ev
GRMATNAME := ${BASNAME}-sl

OUTDIR  := ${BASDIR}/USOL/${RHSFUNC}
OUTNAME := ${OUTDIR}/s

clean-single: 
	-/bin/rm -f ${OUTNAME}.erp  ${OUTNAME}-app.sfn \
          ${OUTNAME}-app-v.eps ${OUTNAME}-app-f.eps \
          ${OUTNAME}-err-v.eps ${OUTNAME}-err-f.eps \
          ${OUTNAME}.ps 
  
run-single: make-approx show-ps

make-approx: ${OUTNAME}-app.sfn

${OUTDIR}:
	mkdir -p ${OUTDIR}

${OUTNAME}-app.sfn: ${PROGDIR}/SPUniSolve Makefile \
            ${OUTDIR} ${BASNAME}.bas \
            ${EVMATNAME}.mat ${GRMATNAME}.mat
	${PROGDIR}/SPUniSolve \
            -rhsName ${RHSFUNC} -coeff ${COEFF} -solName ${SOLNAME} \
            -basisName ${BASNAME} \
            -outName ${OUTNAME} \
            -linearSys cholesky -residual 1 \
            -nonLinearSys iterative 40 1.0e-12 1.0e-12 \
            -smpOrder ${SMPORDER} \
            -plotFinal \
            -paperSize a3 -figSize 80 -meshSize 1.5 \
            -caption "${RUNARGS}" \
          2>&1 | tee ${OUTNAME}.log

# PSVIEW := kghostview
PSVIEW := evince

show-ps: ${OUTNAME}.ps
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-{app,err}-{f,v}.eps` ; do \
          ${PSVIEW} $$el ; \
        done

endif
