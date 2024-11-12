# Last edited on 2023-02-04 07:02:39 by stolfi

SOLVEPROG := SPTimeSpaceSolveG1C0

PROG := ${SOLVEPROG}

PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDERSPACE := 40
SMPORDERTIME := 5

# diffuse:drift:decay:triang:sbasis:fIni:tStep:nSteps:maxIter:RHS:solution

# 0.3333333333333333 4.7e-10
# 0.1666666666666667 2.8e-03
# 0.1631899167899991 3.2e-03
# 0.1546386593239701 4.4e-03
# 0.1538461538461538 4.5e-03
# 0.1428571428571429 6.9e-03
# 0.1540638120433915 

# 0.1540638120433915 0.03125 3.0e-04
# 0.1540638120433915 0.03000 2.3e-04
# 0.1540638120433915 0.02900 3.5e-04
# 0.1540638120433915 0.02800 2.7e-04
# 0.1540638120433915 0.02500 2.3e-05
# 0.1540638120433915 0.02460 1.3e-05
# 0.1540638120433915 0.02445 1.4e-05 (?)
# 0.1540638120433915 0.02435 1.7e-05 (-)
# 0.1540638120433915 0.02400 2.9e-05 (-)
# 0.1540638120433915 0.02000 2.0e-04 (-)

RUNS := \
  +0.1666666666666667:+0.03125:00.00:icosa/raw:nh-d6c1:proc/hrm2:0.100:0300:01:TSFzero:hspn

NONRUNS := \
  +0.1666666666666667:+0.03125:00.00:icosa/raw:nh-d6c1:proc/hrm2:1.000:0005:01:TSFzero:hspn \
  +0.1666666666666667:00.00000:00.00:icosa/raw:nh-d6c1:proc/hrm2:0.100:0005:01:TSFzero:hent \
  +0.1666666666666667:+0.03125:00.00:icosa/raw:nh-d6c1:proc/hrm2:0.100:0005:01:TSFzero:hspn

RHS_RUNS :=

TSTEP_RUNS :=

.PHONY:: run clean-single run-single \
  compute-sol show-ps

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

FIELDS   := ${subst :, ,${RUNARGS}}
DIFFUSE  := ${word  1,${FIELDS}}
DRIFT    := ${word  2,${FIELDS}}
DECAY    := ${word  3,${FIELDS}}
TRIANG   := ${word  4,${FIELDS}}
SBAS     := ${word  5,${FIELDS}}
FINI     := ${word  6,${FIELDS}}
TSTEP    := ${word  7,${FIELDS}}
NSTEPS   := ${word  8,${FIELDS}}
MAXITER  := ${word  9,${FIELDS}}
RHSNAME  := ${word 10,${FIELDS}}
SOLNAME  := ${word 11,${FIELDS}}

TMIN := 0.00000
TMAX := ${shell echo '${TMIN} + ${NSTEPS}*${TSTEP}' | bc -lq}

SBASNAME := ${DATADIR}/sp/pw/${TRIANG}/${SBAS}/b
ININAME  := ${DATADIR}/sp/pw/${TRIANG}/${SBAS}/APP/${FINI}/f-app
OUTDIR   := ${DATADIR}/sp/pw/${TRIANG}/${SBAS}/TSSOL/${RHSNAME}/${FINI}
OUTNAME := ${OUTDIR}/s

run-single: compute-sol show-ps

PSOPTIONS := \
  -ps -paperSize letter \
  -figSize 70 -meshSize 0.75 \
  -caption '%I t = %T' \
  -caption 'step = %D'

EPSOPTIONS := \
  -eps \
  -figSize 45 -meshSize 0.66 \
  -caption '%I t = %T'
  
PLOT_OPTIONS := ${EPSOPTIONS}

compute-sol: ${INIFRAME} ${OUTDIR} ${PROGDIR}/${SOLVEPROG} Makefile
	${PROGDIR}/${SOLVEPROG} \
            -rhsName ${RHSNAME} \
            -cDiffuse '${DIFFUSE}' -cDrift '${DRIFT}' -cDecay '${DECAY}' \
            -iniName ${ININAME} \
            -solName ${SOLNAME} \
            -timeIni '${TMIN}' -timeStep '${TSTEP}' -nSteps ${NSTEPS} \
            -spcBasis ${SBASNAME} \
            -linearSys SVD \
            -nonLinearSys iterative ${MAXITER} '1.0e-12' '1.0e-12' \
            -smpOrder ${SMPORDERSPACE} ${SMPORDERTIME} \
            -verbose \
            -outName ${OUTNAME} \
            ${PLOT_OPTIONS} \
            -plotEvery 10 -writeEvery 10 -printEvery 10

# show-plot -show ${OUTNAME} ${TMIN} ${TMAX}

${OUTDIR}:
	mkdir -p ${OUTDIR}

# PSVIEW := kghostview
PSVIEW := evince

show-ps: ${OUTNAME}.ps
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-{app,err}-{f,v}.eps` ; do \
          ${PSVIEW} $$el ; \
        done

clean-single: 
	-/bin/rm -f ${OUTNAME}.erp \
          ${OUTNAME}.ps 

endif
