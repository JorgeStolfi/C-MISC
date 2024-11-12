# Last edited on 2008-06-10 20:43:18 by stolfi

PROGS := \
  SPComputeMatrices
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 40
SO := ${SMPORDER}

MINVAL := 1.0e-10
MINSV := 0.0

SHOW := 1

# CHECKOPTS := -checkAngles -ignoreTriangs
CHECKOPTS :=

# FBASIS:GBASIS:SMPORDER:EVORTHO

RUNS := \
  pw/icosa/r57/zy/01/h-d3cn-NO:@:${SO}:@ \
  pw/icosa/r57/zy/01/nh-d2cn-NO:@:${SO}:@

RUNS_HARMONIC :=  \
  gl/harm/d004:@:${SO}:-ortho \
  gl/harm/d005:@:${SO}:-ortho \
  gl/harm/d006:@:${SO}:-ortho \
  gl/harm/d012:@:${SO}:-ortho \
  gl/harm/d016:@:${SO}:-ortho \
  gl/harm/d017:@:${SO}:-ortho \
  gl/harm/d018:@:${SO}:-ortho

RUNS_TO_DO := \
  pw/icosa/raw/h-d8c1:@:${SO}:-semi \
  pw/icosa/r57/zy/01/nh-d6c1:@:${SO}:-semi \
   \
  pw/icosa/raw/h-d7c1-NO:@:${SO}:@ \
   \
  pw/icosa/raw/h-d7c1-OLD-NO:@:${SO}:@ \
   \
  pw/icosa/raw/nh-d6c1-NO:@:${SO}:@ \
   \
  pw/icosa/raw/nhx-d6c1-NO:@:${SO}:@ \
   \
  pw/icosa/r57/zz/01/nh-d6c1-NO:@:${SO}:@ \

RUNS_DONE := \
  pw/icosa/r57/zy/01/h-d2cn-NO:@:${SO}:@ \
  pw/icosa/r57/zy/01/nh-d1cn-NO:@:${SO}:@ \
   \
  pw/icosa/raw/h-d8c1:@:${SO}:-semi \
   \
  pw/icosa/r57/zz/01/nh-d6c1:@:${SO}:-semi \
   \
  pw/icosa/raw/nh-d6c1:@:${SO}:-semi \
  pw/icosa/raw/nh-d7c1:@:${SO}:-semi \
   \
  pw/icosa/raw/nhx-d6c1:@:${SO}:-semi \
   \
  pw/icosa/raw/h-d7c1:@:${SO}:-semi \
   \
  pw/icosa/raw/h-d7c1-OLD:@:${SO}:-semi \
   \
  pw/icosa/r57/zz/01/nh-d6c1:@:${SO}:-semi

NONRUNS := \
   \
  pw/icosa/raw/h-d5c1:@:${SO}:-semi \
  pw/icosa/raw/h-d6c1:@:${SO}:-semi \
   \
  pw/icosa/raw/nh-d6c1:@:${SO}:-semi \
  pw/icosa/raw/h-d7c1:@:${SO}:-semi \
   \
  pw/icosa/r3/01/nh-d6c1:@:${SO}:-semi \
  pw/icosa/r3/01/h-d7c1:@:${SO}:-semi 

MUSTREDO := \
  pw/icosa/raw/h-d0cn:@:${SO}:-semi \
  pw/icosa/raw/h-d3cn:@:${SO}:-semi \
   \
  pw/icosa/raw/h-d1c0:@:${SO}:-semi \
  pw/icosa/raw/h-d2c0:@:${SO}:-semi \
  pw/icosa/raw/h-d3c0:@:${SO}:-semi \
  pw/icosa/raw/h-d4c0:@:${SO}:-semi \
  pw/icosa/raw/h-d5c0:@:${SO}:-semi \
   \
  pw/icosa/raw/nh-d2c0:@:${SO}:-semi \
  pw/icosa/raw/nh-d3c0:@:${SO}:-semi \
  pw/icosa/raw/nh-d4c0:@:${SO}:-semi \
  pw/icosa/raw/nh-d5c0:@:${SO}:-semi \
   \
  pw/icosa/raw/h-d5c1:@:${SO}:-semi \
  pw/icosa/raw/h-d6c1:@:${SO}:-semi \
  pw/icosa/raw/h-d7c1:@:${SO}:-semi \
   \
  pw/icosa/raw/nh-d6c1:@:${SO}:-semi \
  pw/icosa/raw/nh-d7c1:@:${SO}:-semi \
   \
  pw/icosa/r57/ah/01/h-d4c0:@:${SO}:-semi \
   \
  pw/icosa/r57/ah/01/nh-d3c0:@:${SO}:-semi \
   \
  pw/icosa/r57/ah/01/h-d7c1:@:${SO}:-semi \
   \
  pw/icosa/r57/ah/01/nh-d6c1:@:${SO}:-semi \
   \
  gl/bez/h-d0:@:${SO}:@ \
  gl/bez/h-d1:@:${SO}:@ \
  gl/bez/h-d2:@:${SO}:@ \
  gl/bez/h-d3:@:${SO}:@ \
  gl/bez/h-d4:@:${SO}:@ \
  gl/bez/h-d5:@:${SO}:@ \
  gl/bez/h-d6:@:${SO}:@ \
  gl/bez/h-d7:@:${SO}:@ \
   \
  gl/bez/nh-d0:@:${SO}:@ \
  gl/bez/nh-d1:@:${SO}:@ \
  gl/bez/nh-d2:@:${SO}:@ \
  gl/bez/nh-d3:@:${SO}:@ \
  gl/bez/nh-d4:@:${SO}:@ \
  gl/bez/nh-d5:@:${SO}:@ \
  gl/bez/nh-d6:@:${SO}:@ \
  gl/bez/nh-d7:@:${SO}:@ \
   \
  gl/poly/nh-d1:@:5

NONRUNS := \
  pw/octa/raw/h-d2c0:@:${SO}:-semi \
  pw/octa/raw/h-d3c0:@:${SO}:-semi \

.PHONY:: run \
  clean-ev clean-ev-single run-ev run-ev-single \
  clean-vg clean-vg-single run-vg run-vg-single \
  clean-sl clean-sl-single run-sl run-sl-single

# all: clean run
# all: run
all: run-ev
# all: run-vg
# all: run-sl
# all: run-vg run-sl

run: run-ev run-vg run-sl
clean:: clean-ev clean-vg clean-sl

clean-ev: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" clean-ev-single ; \
        done
  
run-ev: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" run-ev-single ; \
        done
  
clean-vg: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" clean-vg-single ; \
        done
  
run-vg: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" run-vg-single ; \
        done
  
clean-sl: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" clean-sl-single ; \
        done
  
run-sl: 
	for args in ${RUNS} ; do \
          ${MAKE} RUNARGS="$$args" run-sl-single ; \
        done

  
# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS    := ${subst :, ,${RUNARGS}}
FBASIS    := ${word 1,${FIELDS}}
GBASIS    := ${subst @,${FBASIS},${word 2,${FIELDS}}}
SMPORDER  := ${word 3,${FIELDS}}
EVORTHO   := ${subst @,,${word 4,${FIELDS}}}

FBASNAME := ${DATADIR}/sp/${FBASIS}/b
GBASNAME := ${DATADIR}/sp/${GBASIS}/b

ifeq "/${FBASIS}" "/${GBASIS}"
MATDIR := ${DATADIR}/sp/${FBASIS}
MATNAME := ${MATDIR}/b
BASISOPT := -basis ${FBASNAME}
else
MATDIR := ${DATADIR}/sp/${FBASIS}/x/${GBASIS}
MATNAME := ${MATDIR}/bxb
BASISOPT := -bases ${FBASNAME} ${GBASNAME}
endif

clean-ev-single: 
	-/bin/rm -f ${MATNAME}-ev{,-svdL,-svdR,-choL}.{mat,pgm} ${MATNAME}-ev-svdD.vec
	-/bin/rm -f ${MATNAME}-ev.log
  
clean-vg-single: 
	-/bin/rm -f ${MATNAME}-vg.{mat,pgm}
	-/bin/rm -f ${MATNAME}-vg.log
  
clean-sl-single: 
	-/bin/rm -f ${MATNAME}-sl.{mat,pgm}
	-/bin/rm -f ${MATNAME}-sl.log
  
run-ev-single: ${MATNAME}-ev.mat

run-vg-single: ${MATNAME}-vg.mat

run-sl-single: ${MATNAME}-sl.mat

${MATDIR}:
	mkdir -p ${MATDIR}

${MATNAME}-ev.mat: ${PROGDIR}/SPComputeMatrices Makefile \
            ${MATDIR} ${FBASNAME}.bas ${GBASNAME}.bas
	time ${PROGDIR}/SPComputeMatrices \
            ${BASISOPT} \
            -outName ${MATNAME} \
            -smpOrder ${SMPORDER} -minVal ${MINVAL} \
            -evMatrix ${EVORTHO} -cholesky -SVD \
            ${CHECKOPTS} \
            -verbose \
          2>&1 | tee ${MATNAME}-ev.log
	if [ ${SHOW} ] ; then \
	  ./show-matrix ${MATNAME}-ev ; \
	  ./show-matrix ${MATNAME}-ev-choL ; \
	  ./show-matrix ${MATNAME}-ev-svdL ; \
	  ./show-matrix ${MATNAME}-ev-svdR ; \
	fi

${MATNAME}-vg.mat: ${PROGDIR}/SPComputeMatrices Makefile \
            ${MATDIR} ${FBASNAME}.bas ${GBASNAME}.bas
	time ${PROGDIR}/SPComputeMatrices \
            ${BASISOPT} \
            -outName ${MATNAME} \
            -smpOrder ${SMPORDER} -minVal ${MINVAL} \
            -vgMatrix \
            ${CHECKOPTS} -verbose \
          2>&1 | tee ${MATNAME}-vg.log
	if [ ${SHOW} ] ; then \
	  ./show-matrix ${MATNAME}-vg ; \
	fi

${MATNAME}-sl.mat: ${PROGDIR}/SPComputeMatrices Makefile \
            ${MATDIR} ${FBASNAME}.bas ${GBASNAME}.bas
	time ${PROGDIR}/SPComputeMatrices \
            ${BASISOPT} \
            -outName ${MATNAME} \
            -smpOrder ${SMPORDER} \
            -minVal ${MINVAL} -minSV ${MINSV} \
            -slMatrix \
            ${CHECKOPTS} -verbose \
          2>&1 | tee ${MATNAME}-sl.log
	if [ ${SHOW} ] ; then \
	  ./show-matrix ${MATNAME}-sl ; \
	fi

endif
