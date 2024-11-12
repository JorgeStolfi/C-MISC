# Last edited on 2023-03-31 04:32:49 by stolfi

PROG := tsp_enum
        
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Display options for plot-tour-costs
SHOW_OPTS = -show -timeout 30

# args are NV-TOUR-DISTR-PDIM-NPMAX-SDIM

FEW_CASES := \
  009-1-GEN_GAUSS-007-00000000-010

TRIVIAL_CASES := \
  003-1-GEN_USPHE-002-00000000-000 \
  003-1-GEN_UUNIT-002-00000000-000 \
  003-1-GEN_BIMOD-002-00000000-000 \
  003-1-GEN_GAUSS-002-00000000-000 \

TINY_CASES := \
  004-1-GEN_USPHE-003-00000000-000 \
  004-1-GEN_UUNIT-003-00000000-000 \
  004-1-GEN_BIMOD-003-00000000-000 \
  004-1-GEN_GAUSS-003-00000000-000
  
SMALL_CASES := \
  007-1-GEN_USPHE-006-00000000-000 \
  007-1-GEN_UUNIT-006-00000000-000 \
  007-1-GEN_BIMOD-006-00000000-000 \
  007-1-GEN_GAUSS-006-00000000-000

SMALL_1D_CASES := \
  007-1-GEN_USPHE-001-00000000-000 \
  007-1-GEN_UUNIT-001-00000000-000 \
  007-1-GEN_BIMOD-001-00000000-000 \
  007-1-GEN_GAUSS-001-00000000-000

MINOR_CASES := \
  008-1-GEN_USPHE-007-00000000-000 \
  008-1-GEN_UUNIT-007-00000000-000 \
  008-1-GEN_BIMOD-007-00000000-000 \
  008-1-GEN_GAUSS-007-00000000-000
  
TR_CASES := \
  008-1-GEN_USPHE-007-00000000-000 \
  008-1-GEN_UUNIT-007-00000000-000 \
  008-1-GEN_BIMOD-007-00000000-000 \
  008-1-GEN_GAUSS-007-00000000-000 \
  \
  010-1-GEN_USPHE-009-00100000-000 \
  010-1-GEN_UUNIT-009-00100000-000 \
  010-1-GEN_BIMOD-009-00100000-000 \
  010-1-GEN_GAUSS-009-00100000-000 \
  \
  020-1-GEN_USPHE-021-00100000-000 \
  020-1-GEN_UUNIT-021-00100000-000 \
  020-1-GEN_BIMOD-021-00100000-000 \
  020-1-GEN_GAUSS-021-00100000-000 \
  \
  030-1-GEN_USPHE-031-00100000-000 \
  030-1-GEN_UUNIT-031-00100000-000 \
  030-1-GEN_BIMOD-031-00100000-000 \
  030-1-GEN_GAUSS-031-00100000-000

COST_CASES := \
  009-1-GEN_USPHE-008-00000000-010 \
  009-1-GEN_UUNIT-008-00000000-010 \
  009-1-GEN_BIMOD-008-00000000-010 \
  009-1-GEN_GAUSS-008-00000000-010 \
  009-1-EUC_UCUBE-003-00000000-010 \
  009-1-EUC_UBALL-003-00000000-010 \
  009-1-EUC_USPHE-003-00000000-010 \
  009-1-EUC_GAUSS-003-00000000-010

BIG_CASES := \
  010-1-GEN_UUNIT-002-00100000-000 \
  010-1-GEN_UUNIT-006-00100000-000 \
  010-1-GEN_UUNIT-009-00100000-000 \
  010-1-GEN_UUNIT-015-00100000-000

NON_CASES := \
  \
  004-1-GEN_UUNIT-003-00000000-000 \
  004-1-GEN_UUNIT-003-00000000-000 \
  004-1-GEN_UUNIT-003-00000010-000 \
  004-1-GEN_UUNIT-003-00000010-000 \
  \
  004-1-GEN_UUNIT-002-00000000-000 \
  \
  008-1-GEN_UUNIT-001-00000000-000 \
  008-1-GEN_UUNIT-002-00000000-000 \
  008-1-GEN_UUNIT-007-00000000-000 \
  008-1-GEN_UUNIT-007-00000000-000 \
  008-1-GEN_BIMOD-007-00000000-000 \
  008-1-GEN_BIMOD-007-00000000-000 \
  \
  009-1-GEN_UUNIT-002-00100000-000 \
  009-1-GEN_UUNIT-005-00100000-000 \
  009-1-GEN_UUNIT-008-00100000-000 \
  009-1-GEN_UUNIT-015-00100000-000 

# CASES := ${TR_CASES}
CASES := ${COST_CASES}
# CASES := ${TINY_CASES}
# CASES := ${SMALL_CASES}
# CASES := ${TR_CASES}
# CASES := ${FEW_CASES}

.SUFSDIMS:

.PHONY:: \
  clean-all-dat-did \
  all-cases \
  single-case single-case-both-asymps single-asymp \
  single-cvector

all: clean-all-dat-did all-cases

all-cases:
	for cs in ${CASES}; do \
	  ${MAKE} CASEARGS=$$cs single-case ; \
        done
        
clean-all-dat-did:
	/bin/touch DUMMY.dat DUMMY.did
	/bin/rm -f *.dat *.did

######################################################################
# For recursive make - caller must define ${CASEARGS}

CASEARGS := CASEARGS.IS.UNDEFINED

ifneq "/${CASEARGS}" "/CASEARGS.IS.UNDEFINED"

CASEFIELDS := ${subst -, ,${CASEARGS}}
NV :=        ${word 1,${CASEFIELDS}}
TOUR :=      ${word 2,${CASEFIELDS}}
DISTR :=     ${word 3,${CASEFIELDS}}
PDIM :=      ${word 4,${CASEFIELDS}}
NPMAX :=     ${word 5,${CASEFIELDS}}
SDIM :=      ${word 6,${CASEFIELDS}}

CASENAME := v${NV}-t${TOUR}-${DISTR}-d${PDIM}-m${NPMAX}-f${SDIM}
CVECTORS := 1 2 3 4
  
DAT_DID := out/${CASENAME}.did

${DAT_DID}: ${PROGDIR}/${PROG}
	@echo "=== making .dat files ==="
	for s in ${CVECTORS}; do \
          ${MAKE} CASEARGS=${CASEARGS} CVECTOR=$$s single-cvector ; \
        done
	touch ${DAT_DID}
	@echo "=== done making .dat files ==="

# Global plot with a single cost vector:
GLOB_ONE_EPS := out/${CASENAME}-glob-one.eps

# single-case: ${GLOB_ONE_EPS} single-case-both-asymps
single-case: ${GLOB_ONE_EPS}

${GLOB_ONE_EPS}: ${DAT_DID} plot_all_tour_values.sh
	plot_all_tour_values.sh \
          -nv ${NV} -maxplot 2000 \
          -datastyle linespoints \
          ${SHOW_OPTS} \
          ${CASENAME} ${word 1,${CVECTORS}} \
          > ${GLOB_ONE_EPS}

single-case-both-asymps: 
	for side in lo hi ; do \
          ${MAKE} CASEARGS=${CASEARGS} SIDE=$$side single-asymp ; \
        done
        
######################################################################
# For recursive make - caller must define ${SIDE} as "lo" or "hi"

SIDE := SIDE.IS.UNDEFINED

ifneq "/${SIDE}" "/SIDE.IS.UNDEFINED"

# Simultaneous asymptotic plots with all cost vectors:
RAW_EPS := out/${CASENAME}-asym-${SIDE}-raw.eps
LOG_EPS := out/${CASENAME}-asym-${SIDE}-log.eps
        
single-asymp: ${RAW_EPS} ${LOG_EPS}

# Max number of tours to plot in asymptotic plots
N_ASYMP_PLOT := 200

${RAW_EPS}: ${DAT_DID} plot_small_tour_values.sh
	plot_small_tour_values.sh \
          -nv ${NV} -maxsmall ${N_ASYMP_PLOT} \
          ${SHOW_OPTS} \
          ${CASENAME} ${SIDE} ${CVECTORS} \
          > ${RAW_EPS}

${LOG_EPS}: ${DAT_DID} plot_small_tour_values.sh
	plot_small_tour_values.sh \
          -nv ${NV} -maxsmall ${N_ASYMP_PLOT} \
          -logscale -datastyle linespoints \
          ${SHOW_OPTS} \
          ${CASENAME} ${SIDE} ${CVECTORS} \
          > ${LOG_EPS}

endif
# end ${CVECTOR}
######################################################################
        
######################################################################
# For recursive make - caller must define ${CVECTOR}

CVECTOR := CVECTOR.IS.UNDEFINED

ifneq "/${CVECTOR}" "/CVECTOR.IS.UNDEFINED"

OUTNAME := out/${CASENAME}-c${CVECTOR}

DATFILE := ${OUTNAME}.dat

single-cvector: ${DATFILE}

# Max number of tours to use in asymptotic analysis
N_ASYMP := 100

${DATFILE}: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -seed 4615${CVECTOR}1703 \
            -nv ${NV} -tour ${TOUR} \
            -npmax ${NPMAX} -nsmax ${N_ASYMP} \
            -distr ${DISTR} ${PDIM} \
            -sdim ${SDIM} \
            -outName ${OUTNAME}

endif
# end ${CVECTOR}
######################################################################

endif
# end ${CASEARGS}
######################################################################
