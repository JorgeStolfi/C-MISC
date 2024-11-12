# Last edited on 2023-02-04 07:02:54 by stolfi

PROG := SPComputeApprox

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 40

# FUNCNAME:TRIANG:BASNAME

RUNS := \
  proc/cubx:@:pw/icosa/r57/zy/01/h-d2cn-NO \
  proc/cubx:@:pw/icosa/r57/zy/01/nh-d1cn-NO

RUNS_ANARMHONIC := \
  harm/d012p007:@:pw/icosa/r57/zz/01/nh-d6c1 \
   \
  harm/d005o000:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d005o000:@:pw/icosa/raw/nh-d6c1 \
  harm/d005o000:@:pw/icosa/raw/nhx-d6c1 \
   \
  harm/d005p005:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d005p005:@:pw/icosa/raw/nh-d6c1 \
  harm/d005p005:@:pw/icosa/raw/nhx-d6c1 \
   \
  harm/d012p007:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d012p007:@:pw/icosa/raw/nh-d6c1 \
  harm/d012p007:@:pw/icosa/raw/nhx-d6c1 \
   \
  harm/d006p003:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d006p003:@:pw/icosa/raw/nh-d6c1 \
  harm/d006p003:@:pw/icosa/raw/nhx-d6c1 \
   \
  harm/d006p006:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d006p006:@:pw/icosa/raw/nh-d6c1 \
  harm/d006p006:@:pw/icosa/raw/nhx-d6c1 \
   \
  harm/d007p007:@:pw/icosa/raw/h-d7c1-OLD \
  harm/d007p007:@:pw/icosa/raw/nh-d6c1 \
  harm/d007p007:@:pw/icosa/raw/nhx-d6c1 \
  harm/d007p007:@:pw/icosa/r57/zz/01/nh-d6c1

RUNS_DONE :=  \
  proc/qrtx:@:pw/icosa/raw/nh-d6c1 \
  proc/qrtx:@:pw/icosa/raw/nhx-d6c1 \
   \
  proc/lrry:@:pw/icosa/raw/nh-d6c1 \
  proc/lrry:@:pw/icosa/raw/h-d7c1-OLD \
  proc/lrry:@:pw/icosa/r57/zz/01/nh-d6c1

RUNS_HARMONIC := \
  proc/unit:pw/icosa/raw:gl/harm/d006 \
   \
  proc/hrm2:pw/icosa/raw:gl/harm/d006 \
  proc/hrm2:pw/icosa/raw:gl/harm/d018 \
   \
  proc/lrry:pw/icosa/raw:gl/harm/d006 \
  proc/lrry:pw/icosa/raw:gl/harm/d012 \
  proc/lrry:pw/icosa/raw:gl/harm/d016 \
  proc/lrry:pw/icosa/raw:gl/harm/d017 \
  proc/lrry:pw/icosa/raw:gl/harm/d018 \

NORUNS := \
  proc/spir:@:pw/icosa/r3/01/nh-d6c1 \
  proc/spir:@:pw/icosa/r3/01/h-d7c1 \
  \
  proc/expx:@:pw/icosa/r3/01/nh-d6c1 \
  proc/expx:@:pw/icosa/r3/01/h-d7c1 \
  \
  proc/lrry:@:pw/icosa/r3/01/nh-d6c1 \
  proc/lrry:@:pw/icosa/r3/01/h-d7c1 \
   \
  proc/lrry:@:pw/icosa/raw/h-d7c1  \
  proc/lrry:@:pw/icosa/raw/nh-d6c1 \
  proc/lrry:@:pw/icosa/raw/h-d6c1 \
  \
  proc/hrm2:@:pw/icosa/raw/nh-d6c1 \
  proc/hrm2:@:pw/icosa/raw/h-d7c1 \

NONRUNS := \
  proc/qtrx:@:pw/icosa/raw/h-d4c0 \
   \
  proc/cubx:@:pw/icosa/raw/h-nobary-d3cn \
  proc/cubx:@:pw/icosa/raw/h-d3cn \
  proc/cubx:pw/icosa/raw:gl/bez/h-d3 \
   \
  proc/expx:@:pw/icosa/raw/h-d3cn \
  proc/expx:@:pw/icosa/raw/h-d6c1 \
   \
  proc/expx:@:pw/icosa/raw/h-d2c0 \
  proc/expx:@:pw/icosa/raw/h-d3c0 \
  proc/expx:@:pw/icosa/raw/h-d4c0 \
   \
  proc/expx:@:pw/icosa/raw/nh-d3c0 \
   \
  proc/expx:@:pw/icosa/raw/h-d5c1 \
  proc/expx:@:pw/icosa/raw/h-d7c1 \
  proc/expx:@:pw/icosa/raw/nh-d6c1 \
   \
  proc/expx:@:pw/icosa/r57/ah/01/h-d4c0 \
  proc/expx:@:pw/icosa/r57/ah/01/nh-d3c0 \
   \
  proc/expx:pw/icosa/raw:gl/bez/nh-d4 \
  proc/expx:pw/icosa/raw:gl/bez/h-d3 \
  proc/expx:pw/icosa/raw:gl/bez/h-d4 \
   \
  proc/lrry:@:pw/icosa/raw/h-d4c0 \
  proc/lrry:@:pw/icosa/raw/nh-d3c0 \
   \
  proc/lrry:@:pw/icosa/r57/ah/01/h-d4c0 \
  proc/lrry:@:pw/icosa/r57/ah/01/nh-d3c0 \
   \
  proc/lrry:@:pw/icosa/raw/h-d7c1 \
  proc/lrry:@:pw/icosa/raw/nh-d6c1

.PHONY:: run clean-single run-single

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

FIELDS := ${subst :, ,${RUNARGS}}
FUNC   := ${word 1,${FIELDS}}
TRIANG := ${word 2,${FIELDS}}
BASIS  := ${word 3,${FIELDS}}

FUNNAME := ${DATADIR}/fn/${FUNC}/f
BASNAME := ${DATADIR}/sp/${BASIS}/b

ifeq "/${TRIANG}" "/@" 
  TRINAME := @
  TRIOPTION := 
else
  TRINAME := ${DATADIR}/sp/${TRIANG}/t
  TRIOPTION := -triName ${TRINAME}
endif

MATNAME := ${BASNAME}-ev

# LINMETHOD := cholesky
# FACNAME := ${MATNAME}-choL

LINMETHOD := SVD
FACNAME := ${MATNAME}-svdL

OUTDIR := ${DATADIR}/sp/${BASIS}/APP/${FUNC}
OUTNAME := ${OUTDIR}/f

# PAPEROPTIONS := -paperSize a3 -figSize 80 -meshSize 1.50
PAPEROPTIONS := -paperSize letter -figSize 60 -meshSize 1.00

clean-single: 
	-/bin/rm -f ${OUTNAME}.erp  ${OUTNAME}-app.sfn \
          ${OUTNAME}-app-v.eps ${OUTNAME}-app-f.eps \
          ${OUTNAME}-err-v.eps ${OUTNAME}-err-f.eps \
          ${OUTNAME}.ps
  
run-single: make-approx show-ps

make-approx: ${OUTNAME}-app.sfn

${OUTDIR}:
	mkdir -p ${OUTDIR}

${OUTNAME}-app.sfn: ${PROGDIR}/SPComputeApprox Makefile \
            ${OUTDIR} ${FUNNAME}.sfn \
            ${BASNAME}.bas ${MATNAME}.mat ${FACNAME}.mat
	${PROGDIR}/SPComputeApprox \
            -funcName ${FUNNAME} \
            -basisName ${BASNAME} \
            ${TRIOPTION} \
            -outName ${OUTNAME} \
            -linearSys ${LINMETHOD} \
            -smpOrder ${SMPORDER} -verbose \
            -plot \
            ${PAPEROPTIONS} \
            -caption "${RUNARGS}" \
            -caption "step = %D error = %E" \
          2>&1 | tee ${OUTNAME}.log

# PSVIEW := kghostview
PSVIEW := evince

show-ps: ${OUTNAME}.ps
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-{app,err}-{f,v}.eps` ; do \
          ${PSVIEW} $$el ; \
        done

endif
