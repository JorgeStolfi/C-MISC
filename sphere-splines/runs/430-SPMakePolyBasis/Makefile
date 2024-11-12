# Last edited on 2023-02-04 07:03:07 by stolfi

PROGS := \
  SPMakePolyBasis \
  SPOrthoBasis \
  SPPlotBasis

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

# DEGREE:ORTHO:SHOW

RUNS := \
  004:YES:i0,i1,i2,i3,i4,i5
  
NON_RUNS := \
  004:NO:i0,i1,i2,i3,i4,i5 \
  005:NO:i0,i1,i2,i3,i4,i5 \
  006:NO:i0,i1,i2,i3,i4,i5 \
   \
  012:NO:i0,i1,i2,i3,i4,i5 \
   \
  017:NO:i0,i1,i2,i3,i4,i5,i6,i7,i299,i300,i301,i302,i303,i304,i305,i306 \
  018:NO:i0,i1,i299,i300

.PHONY:: \
  run \
  clean-single run-single \
  make-basis plot-basis show-elems

all: run

clean::
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" clean-single ; \
        done
  
run: 
	for args in ${RUNS} ; do \
          ${MAKE} -f Makefile RUNARGS="$$args" run-single ; \
        done


######################################################################
# For recursive "make" only -- caller must define ${RUNARGS}
# 

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS    := ${subst :, ,${RUNARGS}}
DEGREE    := ${word 1,${FIELDS}}
ORTHO     := ${word 2,${FIELDS}}
SHOW      := ${word 3,${FIELDS}}

# Expand element selection:
COMMA := ,
SHOW := ${subst ${COMMA}, @,@${SHOW}}
SHOW := ${subst @i, -index ,${SHOW}}

run-single: make-basis plot-basis show-elems

TRINAME := ${DATADIR}/sp/pw/octa/raw/t
SMPORDER := 40

ifeq "/${ORTHO}" "/NO"
  OUTDIR := ${DATADIR}/sp/gl/poly/d${DEGREE}-NO
  OUTNAME := ${OUTDIR}/b

  ${OUTNAME}.bas: ${PROGDIR}/SPMakePolyBasis
	mkdir -p ${OUTDIR}
	${PROGDIR}/SPMakePolyBasis \
          -degree ${DEGREE} \
          -outName ${OUTNAME} \
          2>&1 | tee ${OUTNAME}.log
else
  INDIR := ${DATADIR}/sp/gl/poly/d${DEGREE}-NO
  INNAME := ${INDIR}/b

  OUTDIR := ${DATADIR}/sp/gl/poly/d${DEGREE}
  OUTNAME := ${OUTDIR}/b

  # ORTHOOPTS := -gram 1 -eigen 1
  ORTHOOPTS := -gram 1

  ${OUTNAME}.bas: ${PROGDIR}/SPOrthoBasis Makefile \
            ${INNAME}.bas ${OUTDIR} ${TRINAME}.tri
	${PROGDIR}/SPOrthoBasis \
            -basisName ${INNAME} \
            -triName ${TRINAME} \
            ${ORTHOOPTS} -weirdToo \
            -smpOrder ${SMPORDER} \
            -checkOrtho 1.0e-10 \
            -outName ${OUTNAME} \
          2>&1 | tee ${OUTNAME}.log
endif

make-basis: ${OUTNAME}.bas

${OUTDIR}:
	mkdir -p ${OUTDIR}

plot-basis: ${PROGDIR}/SPPlotBasis \
            ${OUTNAME}.bas Makefile
	${PROGDIR}/SPPlotBasis \
          -basisName ${OUTNAME} \
          ${SHOW} \
          -obs 1 10 10 10 \
          -ps -paperSize letter -figSize 40 -meshSize 1.0 \
          -outName ${OUTNAME}

# PSVIEW := kghostview
PSVIEW := evince

show-elems: 
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-el0*.eps` ; do \
          ${PSVIEW} $$el ; \
        done

endif
# end ${RUNARGS}
######################################################################
